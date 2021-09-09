# Statistical Control

# Global parameters and paths
$script_path = "code\GproDIA\src"
$TRIC_path = "$env:ProgramFiles\Anaconda3\Scripts\feature_alignment.py"

<# 
# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$library_name = "fissionyeast"
$library_pqp_file = "$library_name.PQP"
$result_name = $library_name
$result_osw_file = "$result_name.osw"
$result_tsv_file = "$($result_osw_file.Replace('.osw', '')).tsv"
$enable_glycoform_inference = $false
$threads = -1
$test = $true
#>

# Serum Data
$data_path = "data\GproDIA\serum"
$library_name = "serum"
$library_pqp_file = "$($library_name)_uis.PQP"
$result_name = $library_name
$result_osw_file = "$($result_name)_uis.osw"
$result_tsv_file = "$($result_osw_file.Replace('.osw', '')).tsv"
$enable_glycoform_inference = $true
$threads = -1
$test = $true


# Statistical control
cd "$data_path\result\$result_name"

pyprophet merge `
--template="..\..\library\$library_name\$library_pqp_file" `
--out="$result_osw_file" `
(ls *.osw -Name)

python "$script_path\score_glycopeptide_peakgroups.py" `
--level ms2 `
--threads -1 `
--in "$result_osw_file" `
"--$(if ($test) { 'test' } else { 'no-test' })"

if ($enable_glycoform_inference) {
    python "$script_path\score_feature_glycoform.py" `
    --level ms1 `
    --threads -1 `
    --in "$result_osw_file" `
    "--$(if ($test) { 'test' } else { 'no-test' })"

    python "$script_path\score_feature_glycoform.py" `
    --level transition `
    --threads -1 `
    --in "$result_osw_file" `
    "--$(if ($test) { 'test' } else { 'no-test' })"
    
    python "$script_path\infer_glycoforms.py" `
    --in "$result_osw_file"
}

python "$script_path\infer_glycopeptides.py" `
--context global `
--in "$result_osw_file"

if (-not $enable_glycoform_inference) {
    python "$script_path\export_results.py" `
    --in "$result_osw_file" `
    --out "$result_tsv_file" `
    --format legacy_merged `
    --no-glycoform `
    --max_rs_peakgroup_qvalue 0.05 `
    --max_global_glycopeptide_qvalue 0.01 `
    --no-transition_quantification
}
else {
    python "$script_path\export_results.py" `
    --in "$result_osw_file" `
    --out "$result_tsv_file" `
    --format legacy_merged `
    --glycoform --max_glycoform_qvalue 0.05 `
    --max_rs_peakgroup_qvalue 0.05 `
    --max_global_glycopeptide_qvalue 0.01 `
    --no-transition_quantification
}


# TRIC
python $TRIC_path `
--in "$result_tsv_file" `
--out "$($result_tsv_file.Replace('.tsv', ''))_aligned.tsv" `
--file_format openswath `
--fdr_cutoff 0.01 `
--max_fdr_quality 0.2 `
--mst:useRTCorrection True `
--mst:Stdev_multiplier 3.0 `
--method LocalMST `
--max_rt_diff 90 `
--alignment_score 0.001 `
--frac_selected 0 `
--realign_method lowess `
--disable_isotopic_grouping

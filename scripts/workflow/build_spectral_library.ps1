# Spectral Library Building From DDA Results

# Global parameters and paths
$script_path = "code\GproDIA\src"

<# 
# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$dda_main_name = "fissionyeast_6h"
$dda_main_psm_file = "pGlycoDB-GP-FDR-Pro.txt"
$dda_rtcalibrate_name = "fissionyeast_1h"
$dda_rtcalibrate_psm_file = "pGlycoDB-GP-FDR-Pro.txt"
$library_name = "fissionyeast"
$library_pqp_file = "$library_name.PQP"
$swath_window_file = "..\swath_window.tsv"
$enable_glycoform_inference = $false
$rtcalibrate_lowess_frac = 0.667 
$rtcalibrate_lowess_it = 0
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000 
#>

# Serum Data
$data_path = "data\GproDIA\serum"
$dda_main_name = "serum_fraction"
$dda_main_psm_file = "pGlycoDB-GP-FDR-Pro.txt"
$dda_rtcalibrate_name = "serum_1h"
$dda_rtcalibrate_psm_file = "pGlycoDB-GP-FDR-Pro.txt"
$library_name = "serum"
$library_pqp_file = "$($library_name)_uis.PQP"
$swath_window_file = "..\swath_window.tsv"
$background_glycan_file = "..\background_glycan.txt"
$enable_glycoform_inference = $true
$rtcalibrate_lowess_frac = 0.667 
$rtcalibrate_lowess_it = 3
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000
$library_max_background_glycan_number = 50


# Extract assays from pGlyco results
mkdir "$data_path\library\$library_name" -f | Out-Null

python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "$data_path\DDA\$dda_main_name\$dda_main_psm_file" `
--mgf (ls "$data_path\DDA\$dda_main_name\" *.mgf | select -ExpandProperty FullName) `
--out "$data_path\library\$library_name\$dda_main_name.assay.pickle" 

python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "$data_path\DDA\$dda_rtcalibrate_name\$dda_rtcalibrate_psm_file" `
--mgf (ls "$data_path\DDA\$dda_rtcalibrate_name" *.mgf | select -ExpandProperty FullName) `
--out "$data_path\library\$library_name\$dda_rtcalibrate_name.assay.pickle" 

# Build consensus assays
python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\$library_name\$dda_main_name.assay.pickle" `
--out "$data_path\library\$library_name\$($dda_main_name)_consensus.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\$library_name\$dda_rtcalibrate_name.assay.pickle" `
--out "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus.assay.pickle"

# Calibrate retention time
python "$script_path\calibrate_assays_rt.py" `
--multiple_runs `
--smooth lowess --lowess_frac $rtcalibrate_lowess_frac --lowess_it $rtcalibrate_lowess_it `
--in "$data_path\library\$library_name\$($dda_main_name)_consensus.assay.pickle" `
--reference "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus.assay.pickle" `
--out "$data_path\library\$library_name\$($dda_main_name)_consensus_rtcalibrated.assay.pickle" `
--out_anchor "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus_rtanchor.assay.pickle"

# Combine assays and remove redundant assays
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action consensus `
--in "$data_path\library\$library_name\$($dda_main_name)_consensus_rtcalibrated.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\$library_name\$($dda_main_name)_consensus_rtcalibrated_consensus.assay.pickle" `
     "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus.assay.pickle" `
--out "$data_path\library\$library_name\$($library_name)_combined.assay.pickle"

# Filter assays
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\$swath_window_file" `
--min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
--in "$data_path\library\$library_name\$($library_name)_combined.assay.pickle"

if ($enable_glycoform_inference) {
    # Generate glycoform UIS assays
    python "$script_path\generate_glycoform_uis_assays.py" `
    --swath_windows "$data_path\$swath_window_file" `
    --min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
    --background_glycans "$data_path\$background_glycan_file" `
    --max_background_glycan_number $library_max_background_glycan_number `
    --in "$data_path\library\$library_name\$($library_name)_combined_filtered.assay.pickle"
}

# Generate decoy assays
python "$script_path\generate_decoy_assays.py" `
--in "$data_path\library\$library_name\$($library_name)_combined_filtered.assay.pickle"

# Convert assays to PQP
if (-not $enable_glycoform_inference) {
    python "$script_path\convert_assays_to_OpenSWATH_library.py" `
    --disable_glycoform_uis `
    --in "$data_path\library\$library_name\$($library_name)_combined_filtered.assay.pickle" `
    (ls "$data_path\library\$library_name" "$($library_name)_combined_filtered*_decoy.assay.pickle" | select -ExpandProperty FullName) `
    --out "$data_path\library\$library_name\$library_pqp_file"
}
else {
    python "$script_path\convert_assays_to_OpenSWATH_library.py" `
    --enable_glycoform_uis `
    --in "$data_path\library\$library_name\$($library_name)_combined_filtered_uis.assay.pickle" `
    (ls "$data_path\library\$library_name" "$($library_name)_combined_filtered*_decoy.assay.pickle" | select -ExpandProperty FullName) `
    --out "$data_path\library\$library_name\$library_pqp_file"
}

# Convert anchor assays to TraML
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\$swath_window_file" `
--min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
--in "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus_rtanchor.assay.pickle"

python "$script_path\convert_assays_to_traml.py" `
--in "$data_path\library\$library_name\$($dda_rtcalibrate_name)_consensus_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\$library_name\$($library_name)_rtanchor.traML"


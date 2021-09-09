# Spectral Library Building: Semi-Empirically Extended Library

# Global parameters and paths
$script_path = "code\GproDIA\src"

<#
# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$library_empirical_name = "fissionyeast"
$library_empirical_main_name = "fissionyeast"
$library_name = "fissionyeast_extended"
$library_pqp_file = "$library_name.PQP"
$swath_window_file = "..\swath_window.tsv"
$enable_glycoform_inference = $false
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000 
$neighbor_number = 3
$glycopeptide_list = $null
#>

# Serum Data
$data_path = "data\GproDIA\serum"
$library_empirical_name = "serum"
$library_empirical_main_name = "serum"
$library_name = "serum_extend"
$library_pqp_file = "$($library_name)_uis.PQP"
$swath_window_file = "..\swath_window.tsv"
$background_glycan_file = "..\background_glycan.txt"
$enable_glycoform_inference = $true
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000
$library_max_background_glycan_number = 50
$neighbor_number = 3
$glycopeptide_list = "labrepo.glycopeptide.csv"


if ([string]::IsNullOrEmpty($glycopeptide_list)) {
    python "$script_path\generate_semiempirical_assays.py" `
    --action interchange `
    --in "$data_path\library\$library_empirical_name\$($library_empirical_main_name_combined)_combined.assay.pickle" `
    --max_peptide_neighbor_number $neighbor_number `
    --max_glycan_neighbor_number $neighbor_number `
    --min_peptide_occurrence $neighbor_number `
    --min_glycan_occurrence $neighbor_number `
    --glycan_key 'composition' `
    --ignore_glycan_site `
    --out "$data_path\library\$library_name\$($library_empirical_main_name_combined)_semiempirical_k$neighbor_number.assay.pickle"
}
else {
    python "$script_path\generate_semiempirical_assays.py" `
    --action from_list `
    --in "$data_path\library\$library_empirical_name\$($library_empirical_main_name_combined)_combined.assay.pickle" `
    --list "$data_path\library\$library_name\$glycopeptide_list" `
    --glycan_key composition `
    --ignore_glycan_site `
    --max_peptide_neighbor_number $neighbor_number `
    --max_glycan_neighbor_number $neighbor_number `
    --min_peptide_occurrence $neighbor_number `
    --min_glycan_occurrence $neighbor_number `
    --out "$data_path\library\$library_name\$($glycopeptide_list.Replace(".glycopeptide.csv", ""))_$($library_empirical_main_name_combined)_semiempirical_k$neighbor_number.assay.pickle"
}


python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\$library_empirical_name\$($library_empirical_main_name_combined)_combined.assay.pickle" `
     (ls "$data_path\library\$library_name" "*_semiempirical_k$neighbor_number.assay.pickle" | select -ExpandProperty FullName) `
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


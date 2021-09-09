# Spectral Library Building: Mixed Library

# Global parameters and paths
$script_path = "code\GproDIA\src"

# Yeast Serum Mixure Data
$data_path = "data\GproDIA\mix"

$library_1_project_path = "..\fissionyeast"
$library_1_name = "fissionyeast"
$library_1_rtanchor_name = "fissionyeast_1h"
$library_2_project_path = "..\serum"
$library_2_name = "serum"
$library_2_rtanchor_name = "serum_1h"

$library_name = "yeast_serum"
$library_pqp_file = "$($library_name)_uis.PQP"
$swath_window_file = "..\swath_window.tsv"
$background_glycan_file = "..\background_glycan.txt"
$enable_glycoform_inference = $true
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000
$library_max_background_glycan_number = 50



mkdir "$data_path\library\$library_name" -f | Out-Null

if ($enable_glycoform_inference) {
    # Generate glycoform UIS assays
    python "$script_path\generate_glycoform_uis_assays.py" `
    --swath_windows "$data_path\$swath_window_file" `
    --min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
    --background_glycans "$data_path\$background_glycan_file" `
    --max_background_glycan_number $library_max_background_glycan_number `
    --in "$library_1_project_path\library\$library_1_name\$($library_1_name)_combined_filtered.assay.pickle" `
         "$library_2_project_path\library\$library_2_name\$($library_2_name)_combined_filtered.assay.pickle" `
    --out "$data_path\library\$library_name\$($library_name)_filtered_uis.assay.pickle"
}

# Generate decoy assays
python "$script_path\generate_decoy_assays.py" `
--in "$library_1_project_path\library\$library_1_name\$($library_1_name)_combined_filtered.assay.pickle" `
     "$library_2_project_path\library\$library_2_name\$($library_2_name)_combined_filtered.assay.pickle" `
--out_peptide_decoy "$data_path\library\$library_name\$($library_name)_filtered_peptide_decoy.assay.pickle" `
--out_glycan_decoy "$data_path\library\$library_name\$($library_name)_filtered_glycan_decoy.assay.pickle" `
--out_both_decoy "$data_path\library\$library_name\$($library_name)_filtered_both_decoy.assay.pickle"

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
python "$script_path\convert_assays_to_traml.py" `
--in "$library_1_project_path\library\$library_1_name\$($library_1_rtanchor_name)_consensus_rtanchor_filtered.assay.pickle" `
     "$library_2_project_path\library\$library_2_name\$($library_2_rtanchor_name)_consensus_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\$library_name\$($library_name)_rtanchor.traML"


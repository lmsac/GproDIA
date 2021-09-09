# Spectral Library Building: Entrapment Library

# Global parameters and paths
$script_path = "code\GproDIA\src"

# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$organism_label = "fissionyeast"
$library_organism_name = "fissionyeast"

# Serum Data
$entrapment_label = "human"
$library_entrapment_project_path = "..\serum"
$library_entrapment_name = "serum"
$library_entrapment_main_name = "serum_fraction"
$entrapment_size = 500

$library_group_name = "entrapment"
$swath_window_file = "..\swath_window.tsv"
$background_glycan_file = "..\background_glycan.txt"
$enable_glycoform_inference = $true
$library_min_fragment_mz = 200
$library_max_fragment_mz = 2000
$library_max_background_glycan_number = 50


mkdir "$data_path\library\$library_group_name" -f | Out-Null

python -c @"
import sys
sys.path.append(r'$script_path')

from util import load_pickle, save_pickle

assay_file = r'$data_path\$library_entrapment_project_path\library\$library_entrapment_name\$($library_entrapment_main_name)_consensus_rtcalibrated_consensus.assay.pickle'
out_file = r'$data_path\library\$library_group_name\entrapment_$($entrapment_label).assay.pickle'

assays = load_pickle(assay_file)

assays = [
    assay for assay in assays 
    if assay['glycanStruct'].count('F') + assay['glycanStruct'].count('A') > 1
]

save_pickle(assays, out_file)
print('assays saved: {0}, {1} spectra'.format(out_file, len(assays)))
"@


python "$script_path\generate_semiempirical_assays.py" `
--action exchange `
--peptide_assay "$data_path\library\$library_group_name\entrapment_$($entrapment_label).assay.pickle" `
--glycan_assay "$data_path\library\$library_organism_name\$($library_organism_name)_combined.assay.pickle" `
--out "$data_path\library\$library_group_name\entrapment_$($entrapment_label)_peptide_$($organism_label)_glycan.assay.pickle"

python "$script_path\generate_semiempirical_assays.py" `
--action exchange `
--glycan_assay "$data_path\library\$library_group_name\entrapment_$($entrapment_label).assay.pickle" `
--peptide_assay "$data_path\library\$library_organism_name\$($library_organism_name)_combined.assay.pickle" `
--out "$data_path\library\$library_group_name\entrapment_$($organism_label)_peptide_$($entrapment_label)_glycan.assay.pickle"

ls "$data_path\library\$library_group_name" entrapment_*.assay.pickle | % {
    python "$script_path\filter_assays.py" `
    --swath_windows "$data_path\$swath_window_file" `
    --min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
    --in $_.FullName
}


ls "$data_path\library\entrapment" entrapment_*_filtered.assay.pickle | % {
    python -c @"
import sys
sys.path.append(r'$script_path')

import random
from util import load_pickle, save_pickle

n = $entrapment_size
assay_file = r'$($_.FullName)'
assays = load_pickle(assay_file)
assays = random.sample(assays, n)

out_file = assay_file.replace('entrapment_', 'entrapment_%s_' % n)
save_pickle(assays, out_file)
print('assays saved: {0}, {1} spectra'.format(out_file, len(assays)))
"@
}


$library_entrapment_subset_names = @{}
$library_entrapment_subset_names.Add(
    "$entrapment_label", 
    "entrapment_$($entrapment_size)_$($entrapment_label)"
)
$library_entrapment_subset_names.Add(
    "peptide",
    "entrapment_$($entrapment_size)_$($entrapment_label)_peptide_$($organism_label)_glycan"
)
$library_entrapment_subset_names.Add(
    "glycan",
    "entrapment_$($entrapment_size)_$($organism_label)_peptide_$($entrapment_label)_glycan"
)

foreach ($key in $library_entrapment_subset_names.Keys) {
    $library_entrapment_subset_name = $library_entrapment_subset_names[$key]
    $library_name = "$($library_organism_name)_entrapment_$($key)"

    if ($enable_glycoform_inference) {
        # Generate glycoform UIS assays
        python "$script_path\generate_glycoform_uis_assays.py" `
        --swath_windows "$data_path\$swath_window_file" `
        --min_fragment_mz $library_min_fragment_mz --max_fragment_mz $library_max_fragment_mz `
        --background_glycans "$data_path\$background_glycan_file" `
        --max_background_glycan_number $library_max_background_glycan_number `
        --in "$data_path\library\$library_organism_name\$($library_organism_name)_combined_filtered.assay.pickle" `
            "$data_path\library\$library_group_name\$($library_entrapment_subset_name)_filtered.assay.pickle" `
        --out "$data_path\library\$library_group_name\$($library_name)_filtered_uis.assay.pickle"
    }

    # Generate decoy assays
    python "$script_path\generate_decoy_assays.py" `
    --in "$data_path\library\$library_organism_name\$($library_organism_name)_combined_filtered.assay.pickle" `
        "$data_path\library\$library_group_name\entrapment_$($entrapment_size)_$($entrapment_label)_filtered.assay.pickle" `
    --out_peptide_decoy "$data_path\library\$library_group_name\$($library_name)_filtered_peptide_decoy.assay.pickle" `
    --out_glycan_decoy "$data_path\library\$library_group_name\$($library_name)_filtered_glycan_decoy.assay.pickle" `
    --out_both_decoy "$data_path\library\$library_group_name\$($library_name)_filtered_both_decoy.assay.pickle"

    # Convert assays to PQP
    if (-not $enable_glycoform_inference) {
        python "$script_path\convert_assays_to_OpenSWATH_library.py" `
        --disable_glycoform_uis `
        --in "$data_path\library\$library_group_name\$($library_name)_filtered.assay.pickle" `
        (ls "$data_path\library\$library_group_name" "$($library_name)_filtered*_decoy.assay.pickle" | select -ExpandProperty FullName) `
        --out "$data_path\library\$library_group_name\$($library_name).PQP"
    }
    else {
        python "$script_path\convert_assays_to_OpenSWATH_library.py" `
        --enable_glycoform_uis `
        --in "$data_path\library\$library_group_name\$($library_name)_filtered_uis.assay.pickle" `
        (ls "$data_path\library\$library_group_name" "$($library_name)_filtered*_decoy.assay.pickle" | select -ExpandProperty FullName) `
        --out "$data_path\library\$library_group_name\$($library_name)_uis.PQP"
    }
}


# Semi-Empirical Library Cross Validation

# Global parameters and paths
$script_path = "code\GproDIA\src"

# Fission Yeast Data
$data_path = "data\GproDIA\fissionyeast"
$library_empirical_name = "fissionyeast_lab"
$library_empirical_main_name = "fissionyeast_lab"

$library_name = "semiempirical_validation"
$neighbor_number = 3


python "$script_path\generate_semiempirical_assays.py" `
--action cross_validation `
--in "$data_path\library\$library_empirical_name\$($library_empirical_main_name)_combined.assay.pickle" `
--max_peptide_neighbor_number $neighbor_number `
--max_glycan_neighbor_number $neighbor_number `
--min_peptide_occurrence $neighbor_number `
--min_glycan_occurrence $neighbor_number `
--out "$data_path\library\$library_name\$($library_empirical_main_name)_$($library_name)_k$neighbor_number.assay.pickle"


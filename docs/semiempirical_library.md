# GproDIA Tutorial: Using Semi-Empirical Library

## Prerequisites
Before starting this tutorial, users should build a sample-specific spectral library following **GproDIA Tutorial: Glycoform Inference**.

### System Requirements
In this tutorial, more memory space is needed than that in **GproDIA Tutorial: Glycoform Inference**.

GproDIA has been tested on a workstation with Intel Xeon E5-2690 v3 CPU, 32 GB RAM, and Microsoft Windows Server 2016 Version 1607 (OS Build 14393.2430) operating system.

## Tutorial Data
Tutorial data of a serum sample are available at [ProteomeXchange](http://proteomecentral.proteomexchange.org/)/[iProX](https://www.iprox.org/)  with the data set identifier IPX0002792000. 

### Starting Materials
Starting materials and Python binary spectral libraries in **GproDIA Tutorial: Glycoform Inference** are required in this tutorial.
- library\serum\serum_combined.assay.pickle

Glycopeptide lists of interest are optional but recommended. 
- labrepo.glycopeptide.csv (in serum_extend_library_glycopeptide_list.zip)
- Shu2020.glycopeptide.csv (in serum_extend_library_glycopeptide_list.zip)

`labrepo.glycopeptide.csv` contains glycopeptides collected from a repository of serum in **GproDIA Tutorial: Using Repository-Scale Library**. `Shu2020.glycopeptide.csv` contains glycopeptides collected from a publication on N-linked intact glycopeptides in human serum (Shu et al 2020, doi: [10.1074/mcp.RA119.001791](https://doi.org/10.1074/mcp.ra119.001791)).

## Getting Started 
Open Windows PowerShell, set path to the scripts and the project directory (created in **GproDIA Tutorial: Glycoform Inference**) as global parameters.
``` powershell
$script_path = "Path_to_GproDIA\src"
$data_path = "Path_to_project_data"
```

Create a subdirectory `serum_extend` in the `$data_path\library` directory, and move the glycopeptide list files into the subdirectory.
- library\serum_extend\labrepo.glycopeptide.csv
- library\serum_extend\Shu2020.glycopeptide.csv

## Spectral Library Building
### Semi-Empirical Library Generation
Users can generate semi-empirical MS/MS spectra by swapping and combining the peptide and glycan fragment peaks in empirical spectra of different glycopeptides. However, in such libraries, a significant fraction of glycopeptides are actually not present in the samples of interest at a detectable level, which may compromise detection sensitivity.
``` powershell
# Do NOT run this block!
$neighbor_number = 3
python "$script_path\generate_semiempirical_assays.py" `
--action interchange `
--in "$data_path\library\serum\serum_combined.assay.pickle" `
--max_peptide_neighbor_number $neighbor_number `
--max_glycan_neighbor_number $neighbor_number `
--min_peptide_occurrence $neighbor_number `
--min_glycan_occurrence $neighbor_number `
--glycan_key 'composition' `
--ignore_glycan_site `
--out "$data_path\library\serum_extend\serum_semiempirical_k$neighbor_number.assay.pickle"
```

Instead of enumerating all the peptide-glycan combinations in a spectral library, we suggest researchers focus on a subset of glycopeptides of interest for their specific biological questions. In this tutorial, semi-empirical libraries are generated based on the glycopeptide lists.
``` powershell
$neighbor_number = 3
ls "$data_path\library\serum_lab" *.glycopeptide.csv | % {
    python "$script_path\generate_semiempirical_assays.py" `
    --action from_list `
    --in "$data_path\library\serum\serum_combined.assay.pickle" `
    --list "$data_path\library\serum_lab\$_" `
    --glycan_key composition `
    --ignore_glycan_site `
    --max_peptide_neighbor_number $neighbor_number `
    --max_glycan_neighbor_number $neighbor_number `
    --min_peptide_occurrence $neighbor_number `
    --min_glycan_occurrence $neighbor_number `
    --out "$data_path\library\serum_lab\$($_.Replace(".glycopeptide.csv", ""))_serum_semiempirical_k$neighbor_number.assay.pickle"
}
```
A *k*-nearest neighbor strategy is adapted to extract and combine the peptide and glycan fragment peaks from empirical spectra. The number of nearest neighbors (*k*) and be specified by setting `$neighbor_number`.

### Combining Spectral Libraries
Build a consensus spectral library across runs by removing redundant library entries.
```
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\serum\serum_combined.assay.pickle" `
     "$data_path\library\serum_extend\labrepo_serum_semiempirical_k3.assay.pickle" `
     "$data_path\library\serum_extend\Shu2020_serum_semiempirical_k3.assay.pickle" `
--out "$data_path\library\serum_extend\serum_extend_combined.assay.pickle"
```

### Generating Peptide Query Parameters for OpenSWATH
The combined file (`serum_extend_combined.assay.pickle`)  are converted to PQP format, after filtering, identification transition generation for glycoform inference, and decoy generation as described in **GproDIA Tutorial: Glycoform Inference**.

## Downstream Analysis
Conduct targeted data extraction, statistical control, and multi-run alignment, following **GproDIA Tutorial: Glycoform Inference**.

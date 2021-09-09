# GproDIA Tutorial: Using Repository-Scale Library

## Prerequisites
Before starting this tutorial, users should build a sample-specific spectral library following **GproDIA Tutorial: Glycoform Inference**.

## Tutorial Data
Tutorial data of a serum sample are available at [ProteomeXchange](http://proteomecentral.proteomexchange.org/)/[iProX](https://www.iprox.org/)  with the data set identifier IPX0002792000. 

### Starting Materials
Starting materials and Python binary spectral libraries in **GproDIA Tutorial: Glycoform Inference** are required in this tutorial.
- library\serum\serum_1h_consensus.assay.pickle
- library\serum\serum_fraction_consensus_rtcalibrated_consensus.assay.pickle
- library\serum\serum_1h_consensus_rtanchor_filtered.assay.pickle

A organism/tissue-specific MS/MS spectra repository is also required.
- serum_lab_consensus.assay.pickle (in serum_lab_library.pickle.zip)

## Getting Started 
Open Windows PowerShell, set path to the scripts and the project directory (created in **GproDIA Tutorial: Glycoform Inference**) as global parameters.
``` powershell
$script_path = "Path_to_GproDIA\src"
$data_path = "Path_to_project_data"
```

Create a subdirectory `serum_lab` in the `$data_path\library` directory, and move `serum_lab_consensus.assay.pickle` into the subdirectory.
- library\serum_lab\serum_lab_consensus.assay.pickle

## Spectral Library Building
### Retention Time Calibration
Transform retention times to 1 h gradient.
``` powershell
python "$script_path\calibrate_assays_rt.py" `
--multiple_runs `
--smooth lowess --lowess_frac 0.25 --lowess_it 0 `
--in "$data_path\library\serum_lab\serum_lab_consensus.assay.pickle" `
--reference "$data_path\library\serum\serum_1h_consensus.assay.pickle" `
--out "$data_path\library\serum_lab\serum_lab_consensus_rtcalibrated.assay.pickle" `
--out_anchor "$data_path\library\serum_lab\serum_1h_consensus_rtanchor.assay.pickle"
```
Retention time calibration is visualized in a report file (`*_rtcalibrated.pdf`). If retention time is not calibrated properly, you may change the `--lowess_frac` and `--lowess_it` parameters.

### Combining Spectral Libraries
Build a consensus spectral library across runs by removing redundant library entries.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action consensus `
--in "$data_path\library\serum_lab\serum_lab_consensus_rtcalibrated.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\serum\serum_fraction_consensus_rtcalibrated_consensus.assay.pickle" `
     "$data_path\library\serum_lab\serum_lab_consensus.assay.pickle" `
     "$data_path\library\serum\serum_1h_consensus.assay.pickle" `
--out "$data_path\library\serum_lab\serum_lab_combined.assay.pickle"
```

### Combining Retention Time Anchors
Combine the retention time anchors and convert them  to TraML format.
``` powershell
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--in "$data_path\library\serum_lab\serum_1h_consensus_rtanchor.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\serum\serum_1h_consensus_rtanchor_filtered.assay.pickle" `
     "$data_path\library\serum_lab\serum_1h_consensus_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\serum_lab\serum_1h_combined_rtanchor_filtered.assay.pickle"

python "$script_path\convert_assays_to_traml.py" `
--in "$data_path\library\serum_lab\serum_1h_combined_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\serum_lab\serum_lab_rtanchor.traML"
```

### Generating Peptide Query Parameters for OpenSWATH
The combined file (`serum_lab_combined.assay.pickle`)  are converted to PQP format, after filtering, identification transition generation for glycoform inference, and decoy generation as described in **GproDIA Tutorial: Glycoform Inference**.

## Downstream Analysis
Conduct targeted data extraction, statistical control, and multi-run alignment, following **GproDIA Tutorial: Glycoform Inference**.

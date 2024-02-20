# GproDIA Tutorial: Getting Started

## Prerequisites
### System Requirements
GproDIA has been tested on a workstation with Intel Xeon E5-2690 v3 CPU, 16 GB RAM, and Microsoft Windows Server 2016 Version 1607 (OS Build 14393.2430) operating system.

### Software Dependency
The following software are required:
- Python (version 3.5.6 or later, [Anaconda](https://www.anaconda.com/) distribution is recommended)
- [OpenSWATH](http://openswath.org/) (version 2.6.0)
- [PyProphet](https://github.com/PyProphet/pyprophet) (version 2.1.5)
- [msproteomicstools](https://github.com/msproteomicstools/msproteomicstools) (version 0.11.0)
- [pGlyco](http://pfind.ict.ac.cn/software/pGlyco/index.html) (version 2.2.2 or later)
- MSConvert in [ProteoWizard](http://proteowizard.sourceforge.net/)

GproDIA requires the following Python packages integrated in Anaconda:
- numpy (version 1.18.5)
- pandas (version 0.25.3)
- scipy (version 1.4.1)
- scikit-learn (version 0.22.2.post1)

Later versions may be compatible, but have not been tested.

## Tutorial Data
Tutorial data of a fission yeast sample are available at [ProteomeXchange](http://proteomecentral.proteomexchange.org/)/[iProX](https://www.iprox.org/)  with the data set identifier IPX0002792000. 

### Starting Materials
#### LC-MS/MS Raw Data
LC-MS/MS raw data, including 3 DDA runs of 6 h gradient for spectra library building, 1 DDA run of 1 h gradient for retention time calibration, and 4 DIA technical replicates of 1 h gradient:
- fissionyeast_NGP_DDA_6h_1_201210180855.raw
- fissionyeast_NGP_DDA_6h_2_201211010726.raw
- fissionyeast_NGP_DDA_6h_3.raw
- fissionyeast_NGP_DDA_4.raw
- fissionyeast_NGP_DIA_1.raw
- fissionyeast_NGP_DIA_2.raw
- fissionyeast_NGP_DIA_3.raw
- fissionyeast_NGP_DIA_4.raw

#### SWATH Windows File
- swath_window.tsv

A SWATH windows file should be of tab-separated format, including header:
```
start	end
688	712
712	736
736	760
760	784
...
```

#### Protein Sequence Database
The protein sequence database file (FASTA) of *Schizosaccharomyces pombe* downloaded from Swiss-Prot/[UniProt](https://www.uniprot.org/):
- FissionYeast-S.pombe-SP-1808.fasta

### Saved Results
The result files have also been deposited to the ProteomeXchange/iProX repository.

#### DDA Database Search Results
Database search results of the DDA data, and processed MS/MS sptectra (MGF) output by pGlyco:
- pGlycoDB-GP-FDR-Pro.txt (in fissionyeast_DDA6h_pGlyco.zip)
- fissionyeast_NGP_DDA_6h_1_201210180855_HCDFT.mgf
- fissionyeast_NGP_DDA_6h_2_201211010726_HCDFT.mgf
- fissionyeast_NGP_DDA_6h_3_HCDFT.mgf
- pGlycoDB-GP-FDR-Pro.txt (in fissionyeast_DDA1h_pGlyco.zip)
- fissionyeast_NGP_DDA_4_HCDFT.mgf

#### Spectral Library
A SQLite-based spectral library file (PQP) built from the DDA results, and a transition file (TraML) for retention time normalization:
- fissionyeast.PQP
- fissionyeast_rtanchor.traML

#### DIA Results
A SQLite-based OpenSWATH result file (OSW), and exported text files (TSV) after multi-run alignment:
- fissionyeast.osw
- fissionyeast_aligned.tsv

## Getting Started 
Create a directory for the project.
Open Windows PowerShell, set path to the scripts and the project directory as global parameters.
``` powershell
$script_path = "Path_to_GproDIA\src"
$data_path = "Path_to_project_data"
```

## Spectral Library Building
### Database Searching of DDA Data
Perform database searching on the raw DDA data files using pGlyco. See [pGlyco User Guide](http://pfind.ict.ac.cn/software/pGlyco/pGlyco2%20User%20Guide.pdf) for detailed information.

Create a subdirectory named `DDA` in the project directory. Create subdirectories `fissionyeast_6h` and `fissionyeast_1h` in the `$data_path\DDA` directory, and move the result (`pGlycoDB-GP-FDR-Pro.txt`) and MS/MS spectra files (`*.mgf`) of 6 h and 1 h DDA into the corresponding directories.
- DDA\fissionyeast_6h\pGlycoDB-GP-FDR-Pro.txt
- DDA\fissionyeast_6h\fissionyeast_NGP_DDA_6h_*.mgf
- DDA\fissionyeast_1h\pGlycoDB-GP-FDR-Pro.txt
- DDA\fissionyeast_1h\fissionyeast_NGP_DDA_4_HCDFT.mgf

### Spectrum Annotation
Extract and annotate identified MS/MS spectra from the 6 h DDA results.
``` powershell
mkdir "$data_path\library\fissionyeast" -f | Out-Null

python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "$data_path\DDA\fissionyeast_6h\pGlycoDB-GP-FDR-Pro.txt" `
--mgf (ls "$data_path\DDA\fissionyeast_6h\" "fissionyeast_NGP_DDA_6h_*.mgf" | select -ExpandProperty FullName) `
--out "$data_path\library\fissionyeast\fissionyeast_6h.assay.pickle" 
```
A subdirectory named `library` is created in the project directory. The extracted and annotated spectra are imported into a Python binary file (`*.assay.pickle`).

Perform the same operation for the 1 h DDA results.
``` powershell
python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "$data_path\DDA\fissionyeast_1h\pGlycoDB-GP-FDR-Pro.txt" `
--mgf "$data_path\DDA\fissionyeast_1h\fissionyeast_NGP_DDA_4_HCDFT.mgf" `
--out "$data_path\library\fissionyeast\fissionyeast_1h.assay.pickle" 
```

### Building Consensus Spectra
Combine replicate spectra into consensus spectra within each run.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\fissionyeast\fissionyeast_6h.assay.pickle" `
--out "$data_path\library\fissionyeast\fissionyeast_6h_consensus.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\fissionyeast\fissionyeast_1h.assay.pickle" `
--out "$data_path\library\fissionyeast\fissionyeast_1h_consensus.assay.pickle"
```
The consensus spectra are saved in Python binary files (`*_consensus.assay.pickle`).

### Retention Time Calibration
Transform retention times in 6 h gradient to 1 h gradient.
``` powershell
python "$script_path\calibrate_assays_rt.py" `
--multiple_runs `
--smooth lowess --lowess_frac 0.667 --lowess_it 0 `
--in "$data_path\library\fissionyeast\fissionyeast_6h_consensus.assay.pickle" `
--reference "$data_path\library\fissionyeast\fissionyeast_1h_consensus.assay.pickle" `
--out "$data_path\library\fissionyeast\fissionyeast_6h_consensus_rtcalibrated.assay.pickle" `
--out_anchor "$data_path\library\fissionyeast\fissionyeast_1h_consensus_rtanchor.assay.pickle"
```
The calibrated spectral library is saved as a Python binary file (`*_rtcalibrated.assay.pickle`).  Transitions of the anchors are also saved (`*_rtanchor.assay.pickle`), which will be used as retention time references  in the chromatographic extraction step.

Retention time calibration is visualized in a report file (`*_rtcalibrated.pdf`). If retention time is not calibrated properly, you may change the `--lowess_frac` and `--lowess_it` parameters.

### Combining Spectral Libraries
Build a consensus spectral library across runs by removing redundant library entries.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action consensus `
--in "$data_path\library\fissionyeast\fissionyeast_6h_consensus_rtcalibrated.assay.pickle"
```

Append the 1 h spectral library entries to the 6 h spectral library.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\fissionyeast\fissionyeast_6h_consensus_rtcalibrated_consensus.assay.pickle" `
     "$data_path\library\fissionyeast\fissionyeast_1h_consensus.assay.pickle" `
--out "$data_path\library\fissionyeast\fissionyeast_combined.assay.pickle"
```

### Filtering Library Entries
The transitions can then be optimized using a set of rules.
``` powershell
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--in "$data_path\library\fissionyeast\fissionyeast_combined.assay.pickle"
```
The filtered spectral libraries are saved as Python binary files (`*_filtered.assay.pickle`).

If necessary, the rules for transition selection can be modified with a set of parameters. Use `python "$script_path\filter_assays.py" --help` to see descriptions of the parameters.

### Generating Decoys
Generate decoy libraries.
``` powershell
python "$script_path\generate_decoy_assays.py" `
--in "$data_path\library\fissionyeast\fissionyeast_combined_filtered.assay.pickle"
```
The decoys are saved in Python binary files (`*_decoy.assay.pickle`).

### Generating Peptide Query Parameters for OpenSWATH
Combine the target and decoy spectral libraries and convert them to peptide query parameter (PQP) format.
``` powershell
python "$script_path\convert_assays_to_OpenSWATH_library.py" `
--disable_glycoform_uis `
--in "$data_path\library\fissionyeast\fissionyeast_combined_filtered.assay.pickle" `
 (ls "$data_path\library\fissionyeast" fissionyeast_combined_filtered*_decoy.assay.pickle | select -ExpandProperty FullName) `
--out "$data_path\library\fissionyeast\fissionyeast.PQP"
```
This processed spectral library (including decoys) is saved as a PQP file, which is the input for OpenSWATH.

### Saving Retention Time Anchors for OpenSWATH
Convert the retention time anchors to TraML format.
``` powershell
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\..\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--in "$data_path\library\fissionyeast\fissionyeast_1h_consensus_rtanchor.assay.pickle"

python "$script_path\convert_assays_to_traml.py" `
--in "$data_path\library\fissionyeast\fissionyeast_1h_consensus_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\fissionyeast\fissionyeast_rtanchor.traML"
```
This processed anchors are saved in a TraML files, which is the input for OpenSWATH.

## Targeted Data Extraction
### DIA Data Conversion
Convert the raw DIA data ﬁles to mzML format using MSConvert. See [ProteoWizard Documentation for Users](http://proteowizard.sourceforge.net/doc_users.html) for detailed information. 

Create a subdirectory named `mzML` in the project directory.  Move the mzML files into the `$data_path\mzML` directory.
- mzML\fissionyeast_NGP_DIA_*.mzML

### Running OpenSWATH
Run the OpenSWATH analysis workflow. A tutorial for OpenSWATH is in the [OpenSWATH Documentation](http://openswath.org/en/latest/docs/openswath.html). 
``` powershell
mkdir "$data_path\result\fissionyeast" -f | Out-Null
cd "$data_path\mzML"

ls *.mzML -Name | % {
    OpenSwathWorkflow `
    -Library:retentionTimeInterpretation seconds `
    -RTNormalization:alignmentMethod lowess `
    -mz_extraction_window_ms1 10 -mz_extraction_window_ms1_unit ppm `
    -mz_extraction_window 20 -mz_extraction_window_unit ppm `
    -threads 4 `
    -swath_windows_file "..\swath_window.tsv" `
    -tr "..\library\fissionyeast\fissionyeast.PQP" `
    -tr_irt "..\library\fissionyeast\fissionyeast_rtanchor.traML" `
    -in $_ `
    -out_osw "..\result\fissionyeast\$($_.Replace('.mzML', '.osw'))"
}
```
You can set `-threads` the number of threads allowed to be used by OpenSWATH. Detailed descriptions of the `OpenSwathWorkflow` parameters are in the [OpenMS Documentation](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/UTILS_OpenSwathWorkflow.html).

A subdirectory named `result` is created in the project directory. The OpenSWATH results are saved in SQLite format (`*.osw`).

## Statistical Control
Set the working directory to `$data_path\result\fissionyeast`.
``` powershell
cd "$data_path\result\fissionyeast"
```

### Merging
Merged individual OpenSWATH result for each DIA run by PyProphet.
``` powershell
pyprophet merge `
--template="..\..\library\fissionyeast\fissionyeast.PQP" `
--out="fissionyeast.osw" `
(ls *.osw -Name)
```

### Scoring
Conduct semi-supervised learning and error-rate estimation at peak group level.
``` powershell
python "$script_path\score_glycopeptide_peakgroups.py" `
--level ms2 `
--threads -1 `
--in "fissionyeast.osw" `
--test
```
Note that the `--test` option is used to enable test mode with with fixed random seed to get reproducible results. This option should be turned off in practical scenarios.

### Glycopeptide FDR Estimation
Estimated glycopeptide FDR in the global context.
``` powershell
python "$script_path\infer_glycopeptides.py" `
--context global `
--in "fissionyeast.osw"
```
Glycopeptide inference can also be conducted in run-specific and experiment-wide contexts. See [PyProphet Documentation](http://openswath.org/en/latest/docs/pyprophet.html) for detailed information about contexts. 

### Exporting
Export the results to text report.
``` powershell
python "$script_path\export_results.py" `
--in "fissionyeast.osw" `
--out "fissionyeast.tsv" `
--format legacy_merged `
--max_rs_peakgroup_qvalue 0.05 `
--max_global_glycopeptide_qvalue 0.01 `
--no-transition_quantification
```
The results are saved in tab-separated format (`*.tsv`).

## Multi-run Alignment
Set `$TRIC_path` as the path to TRIC script (`feature_alignment.py`), and run TRIC for multi-run alignment. 
``` powershell
$TRIC_path = "C:\Program Files\Anaconda3\Scripts\feature_alignment.py"

python $TRIC_path `
--in "fissionyeast.tsv" `
--out "fissionyeast_aligned.tsv" `
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
```
The parameters are explained in the [TRIC Tutorial](https://github.com/msproteomicstools/msproteomicstools/blob/master/TRIC-README.md).

The aligned results are saved in tab-separated format (`*_aligned.tsv`).

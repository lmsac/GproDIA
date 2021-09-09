# GproDIA Tutorial: Glycoform Inference

## Prerequisites
See **GproDIA Tutorial: Getting Started**.

## Tutorial Data
Tutorial data of a serum sample are available at [ProteomeXchange](http://proteomecentral.proteomexchange.org/)/[iProX](https://www.iprox.org/)  with the data set identifier IPX0002792000. 

### Starting Materials
#### LC-MS/MS Raw Data
LC-MS/MS raw data, including DDA runs with 20 fractions for spectra library building, 1 DDA run of 1 h gradient for retention time calibration, and 3 DIA technical replicates of 1 h gradient:
- 20200922_serum_F1_1_200925174354.raw
- ...
- 20200922_serum_F20_1.raw
- 20200615_serum_NglycP_1hDDA.raw
- 20200615_serum_NGlycP_1hDIA_rep1.raw
- 20200615_serum_NGlycP_1hDIA_rep2.raw
- 20200615_serum_NGlycP_1hDIA_rep3.raw

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
The protein sequence database file (FASTA) of *Homo sapiens* downloaded from Swiss-Prot/[UniProt](https://www.uniprot.org/):
- Human-H.sapiens-SP-1808.fasta

#### Background Glycans File
- background_glycan.txt

### Saved Results
The result files have also been deposited to the ProteomeXchange/iProX repository.

#### DDA Database Search Results
Database search results of the DDA data, and processed MS/MS sptectra (MGF) output by pGlyco:
- pGlycoDB-GP-FDR-Pro.txt (in serum_DDA_fraction_pGlyco.zip)
- 20200922_serum_F1_1_200925174354_HCDFT.mgf
- ...
- 20200922_serum_F20_1_HCDFT.mgf
- pGlycoDB-GP-FDR-Pro.txt (in serum_DDA1h_singlerun_pGlyco.zip)
- 20200615_serum_NglycP_1hDDA_HCDFT.mgf

#### Spectral Library
A SQLite-based spectral library file (PQP) built from the DDA results, and a transition file (TraML) for retention time normalization:
- serum_uis.PQP
- serum_rtanchor.traML

#### DIA Results
A SQLite-based OpenSWATH result file (OSW), and exported text files (TSV) after multi-run alignment:
- serum_uis.osw
- serum_uis_aligned.tsv

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

Create a subdirectory named `DDA` in the project directory. Create subdirectories `serum_fraction` and `serum_1h` in the `$data_path\DDA` directory, and move the result (`pGlycoDB-GP-FDR-Pro.txt`) and MS/MS spectra files (`*.mgf`) of fractionated and single-run DDA into the corresponding directories.
- DDA\serum_fraction\pGlycoDB-GP-FDR-Pro.txt
- DDA\serum_fraction\20200922_serum_F1_*.mgf
- DDA\serum_1h\pGlycoDB-GP-FDR-Pro.txt
- DDA\serum_1h\20200615_serum_NglycP_1hDDA_HCDFT.mgf

### Spectrum Annotation
Extract and annotate identified MS/MS spectra from the 6 h DDA results.
``` powershell
mkdir "$data_path\library\serum" -f | Out-Null

python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "pGlycoDB-GP-FDR-Pro.txt" `
--mgf (ls "$data_path\DDA\serum_fraction" *.mgf | select -ExpandProperty FullName) `
--out "$data_path\library\serum\serum_fraction.assay.pickle" 
```
A subdirectory named `library` is created in the project directory. The extracted and annotated spectra are imported into a Python binary file (`*.assay.pickle`).

Perform the same operation for the 1 h DDA results.
``` powershell
python "$script_path\build_assays_from_pGlyco_mgf.py" `
--clean_glycan_struct `
--clean_glycan_site `
--psm "pGlycoDB-GP-FDR-Pro.txt" `
--mgf "$data_path\DDA\serum_1h\20200615_serum_NglycP_1hDDA_HCDFT.mgf" `
--out "$data_path\library\serum\serum_1h.assay.pickle" 
```

### Building Consensus Spectra
Combine replicate spectra into consensus spectra within each run.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\serum\serum_fraction.assay.pickle"  `
--out "$data_path\library\serum\serum_fraction_consensus.assay.pickle"

python "$script_path\remove_redundant_assays.py" `
--within_run `
--action consensus `
--in "$data_path\library\serum\serum_1h.assay.pickle" `
--out "$data_path\library\serum\serum_1h_consensus.assay.pickle"
```
The consensus spectra are saved in Python binary files (`*_consensus.assay.pickle`).


### Retention Time Calibration
Transform retention times in 6 h gradient to 1 h gradient.
``` powershell
python "$script_path\calibrate_assays_rt.py" `
--multiple_runs `
--smooth lowess --lowess_frac 0.667 --lowess_it 3 `
--in "$data_path\library\serum\serum_fraction_consensus.assay.pickle" `
--reference "$data_path\library\serum\serum_1h_consensus.assay.pickle" `
--out "$data_path\library\serum\serum_fraction_consensus_rtcalibrated.assay.pickle" `
--out_anchor "$data_path\library\serum\serum_1h_consensus_rtanchor.assay.pickle"
```
The calibrated spectral library is saved as a Python binary file (`*_rtcalibrated.assay.pickle`).  Transitions of the anchors are also saved (`*_rtanchor.assay.pickle`), which will be used as retention time references  in the chromatographic extraction step.

Retention time calibration is visualized in a report file (`*_rtcalibrated.pdf`). If retention time is not calibrated properly, you may change the `--lowess_frac` and `--lowess_it` parameters.

### Combining Spectral Libraries
Build a consensus spectral library across runs by removing redundant library entries.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action consensus `
--in "$data_path\library\serum\serum_fraction_consensus_rtcalibrated.assay.pickle"
```

Append the 1 h spectral library entries to the fractionated spectral library.
``` powershell
python "$script_path\remove_redundant_assays.py" `
--across_run `
--action first `
--glycan_key composition `
--ignore_glycan_site `
--in "$data_path\library\serum\serum_fraction_consensus_rtcalibrated_consensus.assay.pickle" `
     "$data_path\library\serum\serum_1h_consensus.assay.pickle" `
--out "$data_path\library\serum\serum_combined.assay.pickle"
```

### Filtering Library Entries
The transitions can then be optimized using a set of rules.
``` powershell
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\..\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--in "$data_path\library\serum\serum_combined.assay.pickle"
```
The filtered spectral libraries are saved as Python binary files (`*_filtered.assay.pickle`).

If necessary, the rules for transition selection can be modified with a set of parameters. Use `python "$script_path\filter_assays.py" --help` to see descriptions of the parameters.

### Library Generation for Glycoform Inference
Generate a spectral library containing identification transitions for glycoform inference.
``` powershell
python "$script_path\generate_glycoform_uis_assays.py" `
--swath_windows "$data_path\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--background_glycans "$data_path\background_glycan.txt" `
--max_background_glycan_number 50 `
--in "$data_path\library\serum\serum_combined_filtered.assay.pickle"
```
The spectral library is saved as Python binary files (`*_uis.assay.pickle`).

### Generating Decoys
Generate decoy libraries using the target library.
Do NOT use the spectral library containing identification transitions for decoy generation.

``` powershell
python "$script_path\generate_decoy_assays.py" `
--in "$data_path\library\serum\serum_combined_filtered.assay.pickle"
```
The decoys are saved in Python binary files (`*_decoy.assay.pickle`).

### Generating Peptide Query Parameters for OpenSWATH
Combine the target and decoy spectral libraries and convert them to peptide query parameter (PQP) format. Be sure that `--enable_glycoform_uis` option is switched on.
``` powershell
python "$script_path\convert_assays_to_OpenSWATH_library.py" `
--enable_glycoform_uis `
--in "$data_path\library\serum\serum_combined_filtered_uis.assay.pickle" `
 (ls "$data_path\library\serum" serum_combined_filtered*_decoy.assay.pickle | select -ExpandProperty FullName) `
--out "$data_path\library\serum\serum_uis.PQP"
```
This processed spectral library (including decoys) is saved as a PQP file, which is the input for OpenSWATH.

### Saving Retention Time Anchors for OpenSWATH
Convert the retention time anchors to TraML format.
``` powershell
python "$script_path\filter_assays.py" `
--swath_windows "$data_path\swath_window.tsv" `
--min_fragment_mz 200 --max_fragment_mz 2000 `
--in "$data_path\library\serum\serum_1h_consensus_rtanchor.assay.pickle"

python "$script_path\convert_assays_to_traml.py" `
--in "$data_path\library\serum\serum_1h_consensus_rtanchor_filtered.assay.pickle" `
--out "$data_path\library\serum\serum_rtanchor.traML"
```
This processed anchors are saved in a TraML files, which is the input for OpenSWATH.

## Targeted Data Extraction
### DIA Data Conversion
Convert the raw DIA data ﬁles to mzML format using MSConvert. See [ProteoWizard Documentation for Users](http://proteowizard.sourceforge.net/doc_users.html) for detailed information. 

Create a subdirectory named `mzML` in the project directory.  Move the mzML files into the `$data_path\mzML` directory.
- mzML\20200615_serum_NGlycP_1hDIA_*.mzML

### Running OpenSWATH
Run the OpenSWATH analysis workflow. A tutorial for OpenSWATH is in the [OpenSWATH Documentation](http://openswath.org/en/latest/docs/openswath.html). 
``` powershell
mkdir "$data_path\result\serum" -f | Out-Null
cd "$data_path\mzML"

ls *.mzML -Name | % {
    OpenSwathWorkflow `
    -Library:retentionTimeInterpretation seconds `
    -RTNormalization:alignmentMethod lowess `
    -RTNormalization:estimateBestPeptides `
    -mz_extraction_window_ms1 10 -mz_extraction_window_ms1_unit ppm `
    -mz_extraction_window 20 -mz_extraction_window_unit ppm `
    -enable_ms1 true -enable_ipf true `
    -threads 4 `
    -swath_windows_file "..\swath_window.tsv" `
    -tr "..\library\serum\serum_uis.PQP" `
    -tr_irt "..\library\serum\serum_rtanchor.traML" `
    -in $_ `
    -out_osw "..\result\serum\$($_.Replace('.mzML', '.osw'))"
}
```
Note that `-enable_ms1` and `-enable_ipf` options  are set to `true` to enable MS1 and transition-level scoring.

You can set `-threads` the number of threads allowed to be used by OpenSWATH. Detailed descriptions of the `OpenSwathWorkflow` parameters are in the [OpenMS Documentation](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/UTILS_OpenSwathWorkflow.html).

A subdirectory named `result` is created in the project directory. The OpenSWATH results are saved in SQLite format (`*.osw`).

## Statistical Control
Set the working directory to `$data_path\result\serum`.
``` powershell
cd "$data_path\result\serum"
```

### Merging
Merged individual OpenSWATH result for each DIA run by PyProphet.
``` powershell
pyprophet merge `
--template="..\..\library\serum\serum_uis.PQP" `
--out="serum_uis.osw" `
(ls *.osw -Name)
```

### Scoring
Conduct semi-supervised learning and error-rate estimation at peak group level.
``` powershell
python "$script_path\score_glycopeptide_peakgroups.py" `
--level ms2 `
--threads -1 `
--in "serum_uis.osw" `
--test
```
Note that the `--test` option is used to enable test mode with with fixed random seed to get reproducible results. This option should be turned off in practical scenarios.

### Glycoform Inference
Conduct the scoring on MS1 and transition-level after the MS2 peak group scoring. 
``` powershell
python "$script_path\score_feature_glycoform.py" `
--level ms1 `
--threads -1 `
--in "serum_uis.osw" `
--test

python "$script_path\score_feature_glycoform.py" `
--level transition `
--threads -1 `
--in "serum_uis.osw" `
--test
```

Apply glycoform inference after scoring.
``` powershell
python "$script_path\infer_glycoforms.py" `
--in "serum_uis.osw"
```

### Glycopeptide FDR Estimation
Estimated glycopeptide FDR in the global context.
``` powershell
python "$script_path\infer_glycopeptides.py" `
--context global `
--in "serum_uis.osw"
```
Glycopeptide inference can also be conducted in run-specific and experiment-wide contexts. See [PyProphet Documentation](http://openswath.org/en/latest/docs/pyprophet.html) for detailed information about contexts. 

### Exporting
Export the results to text report.
``` powershell
python "$script_path\export_results.py" `
--in "serum_uis.osw" `
--out "serum_uis.tsv" `
--format legacy_merged `
--glycoform --max_glycoform_qvalue 0.05 `
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
--in "serum_uis.tsv" `
--out "serum_uis_aligned.tsv" `
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

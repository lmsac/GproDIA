# GproDIA

GproDIA is a pipeline for characterization of intact glycopeptides from DIA data with comprehensive statistical control.

## Dependency
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

## Tutorial
Tutorials are avaliable in the [`docs`](docs) folder.

### Getting Started
[**GproDIA Tutorial: Getting Started**](docs/getting_started.md) describes the analysis workflow for a fission yeast dataset.

### Glycoform Inference
[**GproDIA Tutorial: Glycoform Inference**](docs/glycoform_inference.md) describes a complete analysis workflow including glycoform inference for a serum dataset.

### Repository-Scale Library
[**GproDIA Tutorial: Using Repository-Scale Library**](docs/repository_scale_library.md) describes the analysis workflow for the serum dataset ultilizing a organism-specific repository of MS/MS spectra.

### Semi-Empirical Library
[**GproDIA Tutorial: Using Semi-Empirical Library**](docs/semiempirical_library.md) describes the analysis workflow for the serum dataset with extended library coverage by generating semi-empirical library spectra.

## Publications
Yang, Y., Yan, G., Kong, S., Wu, M., Yang, P., Cao, W., Qiao, L. GproDIA enables data-independent acquisition glycoproteomics with comprehensive statistical control. *Nat Commun* **12**, 6073 (2021). http://doi.org/10.1038/s41467-021-26246-3

## License
GproDIA is distributed under a BSD license. See the LICENSE file for details.

## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn.
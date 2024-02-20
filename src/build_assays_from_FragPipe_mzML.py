import argparse

parser = argparse.ArgumentParser(
    description='Build assays from FragPipe results and mzML files.'
)
parser.add_argument(
    '--psm', nargs='+',
    help='input FragPipe result files'
)
parser.add_argument(
    '--mzml', nargs='+',
    help='input mzML files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)

parser.add_argument(
    '--glycans',
    help='glycan files'
)
parser.add_argument(
    '--fdr', default=0.01, type=float,
    help='total FDR threshold (default: %(default)s)'
)

args = ['--psm', r'D:\FD_7600_Glyco\fragpipe_CID_60min\psm.tsv',
        '--mzml', r'D:\FD_7600_Glyco\mzML\20230324_FD_GlycoPepide_CID_60min_01-20230323_FD_GlycoPepide_CID_60min_01.mzML',
                  r'D:\FD_7600_Glyco\mzML\20230324_FD_GlycoPepide_CID_60min_02-20230323_FD_GlycoPepide_CID_60min_02.mzML',
        '--glycans', r'D:\GlycoDIA\data\background_glycan.txt']

args = parser.parse_args(args)
psm_report_files = args.psm
spectra_files = args.mzml
out_file = args.out
glycan_file = args.glycans
fdr_cutoff = args.fdr

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('psm_report_files', None) is None:
    psm_report_files = list_files(
        path='.',
        pattern='psm\\.tsv$'
    )

if len(psm_report_files) == 0:
    raise ValueError('no psm report files')

# %%
if globals().get('spectra_files', None) is None:
    spectra_files = list_files(
        path='.',
        pattern='\\.mzML$'
    )

if len(spectra_files) == 0:
    raise ValueError('no spectra files')

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(spectra_files[0])[0]
    if len(spectra_files) > 1:
        out_file += '_' + str(len(spectra_files))
    out_file += '.assay.pickle'

# %%
if globals().get('fdr_cutoff', None) is None:
    fdr_cutoff = 0.01

logging.info('use FDR cutoff: ' + str(fdr_cutoff))

# %%
if globals().get('glycan_file', None) is None:
    glycan_file = list_files(
        path='.',
        pattern='glycans\\.txt$'
    )


# %%
import pandas as pd

from util import save_pickle
from spectra.mzmlreader import MzmlReader
from fragpipe import extract_assays_from_spectra

# %%
logging.info('loading psm report(s): ' + '; '.join(psm_report_files))

psm_report = pd.concat(
    (pd.read_table(f) for f in psm_report_files),
    ignore_index=True
)

logging.info('psm report(s) loaded: {0} spectra' \
    .format(len(psm_report)))

# %%
logging.info('loading glycans: ' + glycan_file)

glycan_struct = pd.read_csv(glycan_file, header=None) \
    [0].values.tolist()

logging.info('glycans loaded: {0} glycans' \
            .format(len(glycan_struct)))

# %%
def get_spectra(spectra_file):
    with MzmlReader(spectra_file) as reader:
        while True:
            spec = reader.read_spectrum()
            if spec is None:
                break
            yield spec

# %%
assays = []

for spectra_file in spectra_files:
    logging.info('converting spectra to assays: ' + spectra_file)

    assay_data = extract_assays_from_spectra(
        psm_report=psm_report,
        spectra=get_spectra(spectra_file),
        glycan_struct=glycan_struct,
        glycan_fdr_cutoff=fdr_cutoff,
        return_generator=False
    )
    assays.extend(assay_data)

    logging.info('assays converted: {0}, {1} spectra' \
        .format(spectra_file, len(assay_data)))

logging.info('assays converted: {0} spectra totally' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))


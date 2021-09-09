import argparse

parser = argparse.ArgumentParser(
    description='Build assays from pGlyco results and MGF files.'
)
parser.add_argument(
    '--psm', nargs='+',
    help='input pGlyco result files'
)
parser.add_argument(
    '--mgf', nargs='+',
    help='input MGF files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)
parser.add_argument(
    '--fdr', default=0.01, type=float,
    help='total FDR threshold (default: %(default)s)'
)

clean_glycan_struct_group = parser.add_mutually_exclusive_group(required=False)
clean_glycan_struct_group.add_argument(
    '--clean_glycan_struct', dest='clean_glycan_struct', action='store_true',
    help='remove suspicious glycan structures (default: %(default)s)'
)
clean_glycan_struct_group.add_argument(
    '--no-clean_glycan_struct', dest='clean_glycan_struct', action='store_false',
    help='do not remove suspicious glycan structures (default: True)'
)
parser.set_defaults(clean_glycan_struct=False)

clean_glycan_site_group = parser.add_mutually_exclusive_group(required=False)
clean_glycan_site_group.add_argument(
    '--clean_glycan_site', dest='clean_glycan_site', action='store_true',
    help='remove suspicious glycan sites (default: %(default)s)'
)
clean_glycan_site_group.add_argument(
    '--no-clean_glycan_site', dest='clean_glycan_site', action='store_false',
    help='do not remove suspicious glycan sites (default: True)'
)
parser.set_defaults(clean_glycan_site=False)

args = parser.parse_args()
psm_report_files = args.psm
spectra_files = args.mgf
out_file = args.out
fdr_cutoff = args.fdr
clean_glycan_struct = args.clean_glycan_struct
clean_glycan_site = args.clean_glycan_site

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
        pattern='pGlycoDB-GP-FDR-Pro\\.txt$'
    )
    
if len(psm_report_files) == 0:
    raise ValueError('no psm report files')

# %%
if globals().get('spectra_files', None) is None:
    spectra_files = list_files(
        path='.', 
        pattern='\\.mgf$'
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
if globals().get('clean_glycan_struct', None) is None:
    clean_glycan_struct = False

logging.info('use clean_glycan_struct: ' + str(clean_glycan_struct))

if globals().get('clean_glycan_site', None) is None:
    clean_glycan_site = False

logging.info('use clean_glycan_site: ' + str(clean_glycan_site))

# %%
import pandas as pd

from util import save_pickle
from spectra.mgfreader import MgfReader
from pglyco import extract_assays_from_spectra

# %%
logging.info('loading psm report(s): ' + '; '.join(psm_report_files))

psm_report = pd.concat(
    (pd.read_table(f) for f in psm_report_files),
    ignore_index=True
)

logging.info('psm report(s) loaded: {0} spectra' \
    .format(len(psm_report)))

# %%
def get_spectra(spectra_file):
    with MgfReader(spectra_file) as reader:
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
        total_fdr_cutoff=fdr_cutoff,
        clean_glycan_struct=clean_glycan_struct,
        clean_glycan_site=clean_glycan_site,
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


import argparse

parser = argparse.ArgumentParser(
    description='Generate identification transitions for glycoform inference.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)

parser.add_argument(
    '--swath_windows',
    help='SWATH isolation window file'
)
parser.add_argument(
    '--background_glycans',
    help='background glycan file'
)

parser.add_argument(
    '--max_background_glycan_number', default=20, type=int,
    help='maximum number of background glycans in a SWATH isolation window (default: %(default)s)'
)

ms2_precursor_group = parser.add_mutually_exclusive_group(required=False)
ms2_precursor_group.add_argument(
    '--enable_identification_ms2_precursors', 
    dest='enable_identification_ms2_precursors', action='store_true',
    help='generate identification transitions for MS2-level precursor ions (default: %(default)s)'
)
ms2_precursor_group.add_argument(
    '--disable_identification_ms2_precursors', 
    dest='enable_identification_ms2_precursors', action='store_false',
    help='do not generate identification transitions for MS2-level precursor ions (default: False)'
)
parser.set_defaults(enable_identification_ms2_precursors=True)

frag_filter_group = parser.add_argument_group('fragment filters') 
frag_filter_group.add_argument(
    '--min_fragment_mz', type=float,
    help='lower m/z limit of fragment ions'
)
frag_filter_group.add_argument(
    '--max_fragment_mz', type=float,
    help='upper m/z limit of fragment ions'
)


args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
swath_window_file = args.swath_windows
background_glycan_file = args.background_glycans
max_background_glycan_number = args.max_background_glycan_number
enable_identification_ms2_precursors = args.enable_identification_ms2_precursors

min_fragment_mz = args.min_fragment_mz
max_fragment_mz = args.max_fragment_mz

# %%
import logging

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if globals().get('assay_files', None) is None:
    assay_files = list_files(
        path='.', 
        pattern='\\.assay\\.pickle$'
    )
    
if len(assay_files) == 0:
    raise ValueError('no assay files')
    
# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    out_file += '_uis.assay.pickle'
    
# %%
from util import save_pickle, load_pickle
import pandas as pd

# %%
assays = []
for assay_file in assay_files:
    logging.info('loading assays: ' + assay_file)  
    
    assay_data = load_pickle(assay_file)
    assays.extend(assay_data)
    
    logging.info('assays loaded: {0}, {1} spectra' \
        .format(assay_file, len(assay_data)))

logging.info('assays loaded: {0} spectra totally' \
    .format(len(assays))) 

# %%     
if swath_window_file is not None:
    logging.info('loading SWATH windows: ' + swath_window_file) 
    
    swath_windows = pd.read_csv(swath_window_file, sep='\t')
    
    logging.info('SWATH windows loaded: {0} windows' \
                 .format(len(swath_windows)))
else:
    swath_windows = None
    
# %%
if background_glycan_file is not None:
    logging.info('loading background_glycans: ' + background_glycan_file) 
    
    background_glycans = pd.read_csv(background_glycan_file, header=None) \
        [0].values.tolist()
    
    logging.info('background_glycans loaded: {0} glycans' \
                 .format(len(background_glycans)))
else:
    background_glycans = None
    
# %%
from assay.glycoformassay import GlycoformUisAssayBuilder

builder = GlycoformUisAssayBuilder(
    background_glycans=background_glycans,
    max_background_glycan_number=max_background_glycan_number,
    enable_identification_ms2_precursors=enable_identification_ms2_precursors
)

logging.info('generating identification transitions') 

assays_uis = builder.build_uis_assays(
    assays, 
    swath_windows=swath_windows,
    min_fragment_mz=min_fragment_mz,
    max_fragment_mz=max_fragment_mz
)

logging.info('identification transitions generated: {0} spectra' \
             .format(len(assays_uis)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays_uis, out_file)

logging.info('assays saved: {0}, {1} spectra' \
             .format(out_file, len(assays_uis)))


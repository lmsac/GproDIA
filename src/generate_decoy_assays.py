import argparse

parser = argparse.ArgumentParser(
    description='Generate decoy assays.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out_peptide_decoy',
    help='output peptide decoy assay file'
)
parser.add_argument(
    '--out_glycan_decoy',
    help='output glycan decoy assay file'
)
parser.add_argument(
    '--out_both_decoy',
    help='output both decoy assay file'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
peptide_decoy_out_file = args.out_peptide_decoy
glycan_decoy_out_file = args.out_glycan_decoy
both_decoy_out_file = args.out_both_decoy
    
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
        pattern='(?<!decoy)\\.assay\\.pickle$'
    )
    
if len(assay_files) == 0:
    raise ValueError('no assay files')
            
# %%
import os

if globals().get('peptide_decoy_out_file', None) is None and \
    globals().get('glycan_decoy_out_file', None) is None and \
    globals().get('both_decoy_out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    peptide_decoy_out_file = out_file + '_peptide_decoy.assay.pickle'
    glycan_decoy_out_file = out_file + '_glycan_decoy.assay.pickle'
    both_decoy_out_file = out_file + '_both_decoy.assay.pickle'
      
# %%
from util import save_pickle, load_pickle
from decoy import GlycoDecoyAssayGenerator

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
decoy = GlycoDecoyAssayGenerator()

# %%
if globals().get('peptide_decoy_out_file', None) is not None:
    logging.info('generating peptide decoy assays')
    
    peptide_decoy_assays = list(map(decoy.peptide_decoy, assays))
    
    logging.info('peptide decoy assays generated: {0} spectra' \
        .format(len(peptide_decoy_assays)))

# %%
#for assay in peptide_decoy_assays:
#    fragments = assay['fragments']
#    if sum([
#        not a.startswith(fragments['fragmentType'][i])
#        for i, a in enumerate(fragments['fragmentAnnotation'])
#    ]) > 0:
#        import pandas as pd
#        print(pd.DataFrame.from_dict(fragments))
#        raise Exception    
# %%   
if 'peptide_decoy_assays' in globals():    
    logging.info('saving peptide decoy assays: {0}' \
        .format(peptide_decoy_out_file))
    
    save_pickle(peptide_decoy_assays, peptide_decoy_out_file)
        
    logging.info('peptide decoy assays saved: {0}, {1} spectra' \
        .format(peptide_decoy_out_file, len(peptide_decoy_assays)))

# %%
if globals().get('glycan_decoy_out_file', None) is not None or \
    globals().get('both_decoy_out_file', None) is not None:
    logging.info('generating glycan decoy assays')
    
    glycan_decoy_assays = list(map(decoy.glycan_decoy, assays))
    
    logging.info('glycan decoy assays generated: {0} spectra' \
        .format(len(glycan_decoy_assays)))

# %%   
if 'glycan_decoy_assays' in globals():    
    logging.info('saving glycan decoy assays: {0}' \
        .format(glycan_decoy_out_file))
    
    save_pickle(glycan_decoy_assays, glycan_decoy_out_file)
        
    logging.info('glycan decoy assays saved: {0}, {1} spectra' \
        .format(glycan_decoy_out_file, len(glycan_decoy_assays)))

# %%
if globals().get('both_decoy_out_file', None) is not None: 
    logging.info('generating both decoy assays')
    
    both_decoy_assays = list(map(decoy.peptide_decoy, glycan_decoy_assays))
    
    logging.info('both decoy assays generated: {0} spectra' \
        .format(len(both_decoy_assays)))

# %%
#for assay in both_decoy_assays:
#    fragments = assay['fragments']
#    if sum([
#        not a.startswith(fragments['fragmentType'][i])
#        for i, a in enumerate(fragments['fragmentAnnotation'])
#    ]) > 0:
#        import pandas as pd
#        print(pd.DataFrame.from_dict(fragments))
#        raise Exception  
    
# %%   
if 'both_decoy_assays' in globals():
    logging.info('saving both decoy assays: {0}' \
        .format(both_decoy_out_file))
    
    save_pickle(both_decoy_assays, both_decoy_out_file)
        
    logging.info('both decoy assays saved: {0}, {1} spectra' \
        .format(both_decoy_out_file, len(both_decoy_assays)))


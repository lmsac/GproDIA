import argparse

parser = argparse.ArgumentParser(
    description='Convert assays to TraML.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output TraML file'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
    
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
    out_file += '.traML'
    
# %%
from util import load_pickle
from openswath.tramlwriter import TramlWriter
from openswath.glycotraml import traml_writer_glyco_parameters
 
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
traml = TramlWriter(parameters=traml_writer_glyco_parameters())

# %%
logging.info('converting assays to TraML')
    
data = traml.assays_to_traml(assays)

logging.info('assays converted: {0} transition groups, {1} transitions' \
    .format(len(traml.peptide_dict), len(traml.transition_dict)))

# %%
logging.info('saving TraML: {0}' \
    .format(out_file))

traml.save_traml(out_file)

logging.info('TraML saved: {0}, {1} transition groups, {2} transitions' \
    .format(out_file, len(traml.peptide_dict), len(traml.transition_dict)))
    
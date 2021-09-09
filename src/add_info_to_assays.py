import argparse

parser = argparse.ArgumentParser(
    description='Add information to assays.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--info', nargs='+',
    help='input information files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)
parser.add_argument(
    '--info_id', nargs='+', default=['file', 'scan'],
    help='information names'
)
parser.add_argument(
    '--info_name', nargs='+',
    help='information names'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
info_files = args.info
out_file = args.out
info_id = args.info_id
info_names = args.info_name

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
if globals().get('info_files', None) is None:
    info_files = list_files(
        path='.', 
        pattern='\\.scaninfo\\.csv$'
    )
    
if len(info_files) == 0:
    raise ValueError('no information files')    

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    out_file += '_' + str(len(assay_files))
    out_file += '.assay.pickle'
    
# %%
from util import load_pickle, save_pickle
import pandas as pd
import re

# %%
logging.info('loading information file(s): ' + '; '.join(info_files))

info_data = pd.concat(
    (
        pd.read_csv(f) if f.endswith('.csv') else pd.read_table(f)
        for f in info_files
    ),
    ignore_index=True
)

logging.info('information file(s) loaded: {0} scans' \
    .format(len(info_data)))

# %%
if globals().get('info_id', None) is None or len(info_id) == 0:
    raise ValueError('no information ID')
   
logging.info('use information id: ' + \
    ', '.join(info_id))

if globals().get('info_names', None) is None:  
    info_names = list(set(info_data.columns) - set(info_id))
    
logging.info('use information: ' + \
    ', '.join(info_names))

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
def get_info_from_assay(assay, name):
    def get_value_from_dict(d, path):
        value = None
        for name in path:
            if isinstance(d, dict):                
                value = d.get(name, None)
                d = value
            else:
                return None                
        return value
                
    value = get_value_from_dict(assay, [name])
    if value is None:
        value = get_value_from_dict(assay, ['metadata', name])    
    if value is None:
        name = name.split('.')
        value = get_value_from_dict(assay, name)
    if value is None:
        value = get_value_from_dict(assay, ['metadata'] + name)

    return value

def set_info_to_assay(assay, name, value):
    def set_value_to_dict(d, path, value):        
        dd = d
        for i, name in enumerate(path):
            if i == len(path) - 1:
                dd[name] = value
            else:
                x = dd.get(name, None)
                if not isinstance(x, dict):
                    x = {}
                    dd[name] = x
                dd = x
        return d
                
    name = name.split('.')
    set_value_to_dict(assay, name, value)    

    return assay

def to_native(x):    
    if re.search('float', type(x).__name__):
        return float(x)
    if re.search('int', type(x).__name__):
        return int(x)
    if re.search('bool', type(x).__name__):
        return bool(x)
    return x
    
# %%
info_id_1 = [
    re.sub('=[^=]+$', '', x)
    for x in info_id
]
info_id_2 = [
    re.sub('^[^=]+=', '', x)
    for x in info_id
]
info_data.set_index(info_id_2, inplace=True)

info_names_1 = [
    re.sub('=[^=]+$', '', x) if '=' in x else 'metadata.' + x
    for x in info_names
]
info_names_2 = [
    re.sub('^[^=]+=', '', x)
    for x in info_names
]

for assay in assays:
    assay_id = tuple(
        get_info_from_assay(assay, x)
        for x in info_id_1
    )   
    info = info_data.loc[assay_id]
    
    for name1, name2 in zip(info_names_1, info_names_2):
         set_info_to_assay(assay, name1, to_native(info[name2]))
    
# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)
    
logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

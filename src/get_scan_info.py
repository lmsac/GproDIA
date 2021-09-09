import argparse

parser = argparse.ArgumentParser(
    description='Get scan information from mzML.'
)

parser.add_argument(
    '--in', 
    help='input mzML file'
)
parser.add_argument(
    '--out',
    help='output CSV file'
)
parser.add_argument(
    '--start_scan', type=int,
    help='start scan number'
)
parser.add_argument(
    '--end_scan', type=int,
    help='end scan number'
)
parser.add_argument(
    '--info_name', nargs='+', 
    default=[
        'msLevel', 'rt', 
        'scanWindowLowerLimit', 'scanWindowUpperLimit',
        'isolationWindowTargetMZ', 
        'isolationWindowLowerOffset', 'isolationWindowUpperOffset',
        'collisionEnergy',
        'compensationVoltage'
    ],
    help='information names'
)

args = parser.parse_args()
mzml_file = getattr(args, 'in')
out_file = args.out
start_scan = args.start_scan
end_scan = args.end_scan
info_names = args.info_name

# %%
import logging

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
if mzml_file is None:
    raise ValueError('no mzML file')
 
# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(mzml_file)[0]
    out_file += '.scaninfo.csv'
    
# %%
if globals().get('info_names', None) is None or len(info_names) == 0:
    raise ValueError('no information names')
      
logging.info('use information: ' + \
    ', '.join(info_names))

# %%
from spectra.mzmlreader import MzmlReader
import csv
import re

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

# %%
info_names_1 = [
    re.sub('=[^=]+$', '', x)
    for x in info_names
]
info_names_2 = [
    re.sub('^[^=]+=', '', x)
    for x in info_names
] 

# %%
reader = MzmlReader(mzml_file)

csv_file = open(out_file, 'w', newline='')
writer = csv.writer(csv_file)

# %%
writer.writerow(['file', 'scan'] + info_names_1)

# %%
logging.info('loading scan information from mzML: ' + mzml_file)
    
scan = 0
while True:
    spec = reader.read_spectrum()    
    if spec is None:
        break
    
    scan += 1
    if start_scan is not None and scan < start_scan:
        continue
    elif end_scan is not None and scan > end_scan:
        break
    
    file = spec['metadata'].get('file', None)
    if file is None:
        title = spec['metadata'].get('title', None)
        mat = re.search('^(.+)\\.([0-9]+)\\.[0-9]+\\.[0-9]*( |$)', title)
        if mat is not None:
            file = mat.group(1)
    if file is None:
        file = os.path.splitext(mzml_file)[0]
        
    scan_number = spec['metadata'].get('scan', None)  
    if scan_number is None:
        title = spec['metadata'].get('title', None)
        mat = re.search('^(.+)\\.([0-9]+)\\.[0-9]+\\.[0-9]*( |$)', title)
        if mat is not None:
            scan_number = mat.group(2)
    if scan_number is None:
        scan_id = spec['metadata'].get('id', None) 
        mat = re.search('scan=([0-9]+)', scan_id)
        if mat is not None:
            scan_number = mat.group(1)
    if scan_number is None:
        scan_number = scan
    scan_number = int(scan_number)
        
    writer.writerow([file, scan_number] + [
        get_info_from_assay(spec, x)
        for x in info_names_2
    ])
    
    if (scan - (start_scan or 1) + 1) % 1000 == 0 or scan == end_scan:
        logging.info(
            'scan {0} processed' \
                .format(scan)
        )
    
# %%
csv_file.close()

logging.info('scan information saved: ' + out_file)

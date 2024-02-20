import argparse

parser = argparse.ArgumentParser(
    description='Update assays from table files.'
)
parser.add_argument(
    '--assay', required=True,
    help='input assay files'
)
parser.add_argument(
    '--data', required=True,
    help='input table files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)

args = parser.parse_args()
assay_file = args.assay
data_file = args.data
out_file = args.out

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)


# %%
import os

if globals().get('out_file', None) is None:
    out_file = assay_file

# %%
import pandas as pd


logging.info('loading data: ' + data_file)

data = pd.read_csv(data_file, sep=',' if data_file.endswith('.csv') else '\t')

logging.info('data loaded: {0} entries' \
    .format(len(data)))


# %%
from util import load_pickle

logging.info('load ions: ' + assay_file)

assays = load_pickle(assay_file)

logging.info('assays loaded: {0}, {1} entries' \
                .format(assay_file, len(assays)))

# %%
from assay.values import assign_assays_values




logging.info('assigning data')
assign_assays_values(
    assays, data=data,
    keys=[{'name': 'sequence', 'path': 'peptideSequence'}],
    params=[{'name': 'protein', 'path': ['metadata', 'protein']}]
)



# %%
from util import save_pickle

logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))


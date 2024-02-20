import argparse

parser = argparse.ArgumentParser(
    description='Convert assays to OpenSWATH spectral library.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output spectral library file'
)
parser.add_argument(
    '--format', choices=['PQP', 'tsv'], default='PQP',
    help='output spectral library format (default: %(default)s)'
)

glycoform_group = parser.add_mutually_exclusive_group(required=False)
glycoform_group.add_argument(
    '--enable_glycoform_uis',
    dest='enable_glycoform_uis', action='store_true',
    help='generate identification transitions for glycoform inference (default: %(default)s)'
)
glycoform_group.add_argument(
    '--disable_glycoform_uis',
    dest='enable_glycoform_uis', action='store_false',
    help='do not generate identification transitions for glycoform inference (default: True)'
)
parser.set_defaults(enable_glycoform_uis=False)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
out_format = args.format
enable_glycoform_uis = args.enable_glycoform_uis

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
if globals().get('out_format', None) is None:
    out_format = 'PQP'

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    out_file += '.' + out_format

# %%
from util import load_pickle
from assay.assay2table import AssayToDataFrameConverter
from openswath import OpenSWATH_glyco_columns

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
assay_to_table = AssayToDataFrameConverter(
    columns=OpenSWATH_glyco_columns(enable_glycoform_uis=enable_glycoform_uis)
)

# %%
logging.info('converting assays to table')

data = assay_to_table.assays_to_dataframe(assays)
data['ProteinId'] = data['ProteinId'].fillna('NA')

logging.info('assays converted: {0} transition groups, {1} transitions' \
    .format(data['TransitionGroupId'].nunique(), len(data)))

# %%
if out_format == 'PQP':
    logging.info('saving PQP: {0}' \
        .format(out_file))

    if enable_glycoform_uis:
        from openswath.glycoformpqp import GlycoformUisQueryParameter
        pqp = GlycoformUisQueryParameter(out_file)
    else:
        from openswath import GlycoPeptideQueryParameter
        pqp = GlycoPeptideQueryParameter(out_file)

    pqp.create_table()
    pqp.add_data(data)
    pqp.close()

    logging.info('PQP saved: {0}, {1} transition groups, {2} transitions' \
        .format(out_file, data['TransitionGroupId'].nunique(), len(data)))

# %%
if out_format == 'tsv':
    logging.info('saving table: {0}' \
        .format(out_file))

    data.to_csv(
        out_file,
        index=False, sep='\t'
    )

    logging.info('table saved: {0}, {1} transition groups, {2} transitions' \
        .format(out_file, data['TransitionGroupId'].nunique(), len(data)))

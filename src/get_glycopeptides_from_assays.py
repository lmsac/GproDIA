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
    out_file += '.glycopeptides.csv'

# %%
from util import load_pickle
from assay.assay2table import AssayToDataFrameConverter
from assay.modseq import ModifiedSequenceConverter

def glycopeptide_columns():
    seq_converter = ModifiedSequenceConverter()

    def modified_sequence(assay, **kwargs):
        r = seq_converter.to_tpp_format(
            sequence=assay['peptideSequence'],
            modification=assay.get('modification', None)
        )
        r = r.replace("C[160]", "cC")
        r = r.replace("M[147]", "oxM")
        r = r.replace("n[43]", "ac")
        r = r.replace("J", "N")
        return r
    columns = [
        {
            'name': 'modified_sequence',
            'function': modified_sequence
        },
        {
            'name': 'glycan_struct',
            'path': 'glycanStruct',
        },
        {
            'name': 'glycan_position',
            'path': 'glycanSite',
        },
        {
            'name': 'precursor_charge',
            'path': 'precursorCharge',
        }
    ]

    columns.extend([
        {
            'name': 'retention_time',
            'path': 'rt'
        },
        {
            'name': 'protein',
            'path': ['metadata', 'protein']
        },
        {
            'name': 'protein_site',
            'path': ['metadata', 'proteinSite']
        },
    ])
    return columns

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
    columns=glycopeptide_columns()
)

# %%
logging.info('converting assays to table')

data = assay_to_table.assays_to_dataframe(assays)

logging.info('assays converted: {0} glycopeptide precursors' \
    .format(len(data)))

# %%
logging.info('saving table: {0}' \
    .format(out_file))

data.to_csv(out_file, index=False)

logging.info('table saved: {0}, {1} glycopeptide precursors' \
    .format(out_file, len(data)))
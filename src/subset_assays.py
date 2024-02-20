import argparse

parser = argparse.ArgumentParser(
    description='Subset assays using a DIA search result.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--result', nargs='+',
    help='input DIA search result files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)


assay_filter_group = parser.add_argument_group('result filters')
assay_filter_group.add_argument(
    '--max_rs_peakgroup_qvalue', type=float, default=0.05,
    help='filter results to maximum run-specific peak group-level q-value (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--max_global_glycopeptide_qvalue', type=float, default=1,
    help='Filter results to maximum global glycopeptide-level q-value (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--max_glycoform_qvalue', type=float, default=0.05,
    help='filter results to maximum glycoform q-value (default: %(default)s)'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
result_files = args.result
out_file = args.out

filter_args = vars(args)
filter_args.pop('in')
filter_args.pop('out')
filter_args.pop('result')

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
if globals().get('result_files', None) is None:
    result_files = list_files(
        path='.',
        pattern='\\.tsv$'
    )

if len(result_files) == 0:
    raise ValueError('no result files')

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    out_file += '_subset.assay.pickle'

# %%
from util import save_pickle, load_pickle
import numpy as np
import pandas as pd

# %%
result = []
for result_file in result_files:
    logging.info('loading result: ' + result_file)

    result_ = pd.read_csv(result_file, sep='\t')
    result.append(result_)

    logging.info('result loaded: {0} records' \
                    .format(len(result_)))

result = pd.concat(result, copy=False, ignore_index=True)

logging.info('result loaded: {0} records totally' \
    .format(len(result)))

# %%
def filter_results(result, max_rs_peakgroup_qvalue=0.05, max_global_glycopeptide_qvalue=0.01, max_glycoform_qvalue=0.05):
    remove = np.zeros(len(result), dtype=bool)
    if 'decoy' in result.columns:
        remove |= result['decoy'] > 0
    if 'm_score' in result.columns and 'ms2_m_score' in result.columns:
        remove |= result['ms2_m_score'] > max_rs_peakgroup_qvalue
        remove |= result['m_score'] > max_glycoform_qvalue
    elif 'm_score' in result.columns:
        remove |= result['m_score'] > max_rs_peakgroup_qvalue
    if 'm_score_glycopeptide_global' in result.columns:
        remove |= result['m_score_glycopeptide_global'] > max_global_glycopeptide_qvalue
    return result.loc[~remove]


logging.info(
    'filtering results using the following parameters: \n' + \
    '\n'.join((
        k + '=' + str(v)
        for k, v in filter_args.items()
        if v is not None
    ))
)
result = filter_results(result, **filter_args)

glycopeptides = result[["FullPeptideName", "GlycanStruct", "GlycanSite", "Charge"]] \
    .drop_duplicates()

logging.info('results filtered: {0} record of {1} glycopeptide precursors remaining' \
    .format(len(result), len(glycopeptides)))


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
from assay.assay2table import AssayToDataFrameConverter
from openswath import OpenSWATH_glyco_columns

assay_to_table = AssayToDataFrameConverter(
    columns=[
        x for x in OpenSWATH_glyco_columns(enable_glycoform_uis=False)
        if x["name"] in {'ModifiedPeptideSequence', 'GlycanStruct', 'GlycanSite', 'PrecursorCharge'}
    ]
)
data = assay_to_table.assays_to_dataframe(assays)
data = data.rename({'ModifiedPeptideSequence': "FullPeptideName", 'PrecursorCharge': "Charge"}, axis=1)
assert len(data) == len(assays)

# %%
data = data.merge(glycopeptides, how="left", indicator=True, copy=False)
assert len(data) == len(assays)

assays = [assay for keep, assay in zip(data["_merge"] == "both", assays) if keep]

logging.info('assays subset: {0} spectra remaining' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

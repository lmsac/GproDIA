import argparse

parser = argparse.ArgumentParser(
    description='Score replicate assays.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--out',
    help='output score file'
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
    out_file += '.score.csv'

# %%
from util import load_pickle

# %%
if assay_files is not None:
    assays = []
    
    for assay_file in assay_files:
        logging.info('loading assays: ' + assay_file)  
        
        assay_data = load_pickle(assay_file)
        assays.extend(assay_data)
        
        logging.info('assays loaded: {0}, {1} spectra' \
            .format(assay_file, len(assay_data)))
    
    logging.info('assays loaded: {0} spectra totally' \
        .format(len(assays))) 
else:
    assays = None

# %%
import pandas as pd
import itertools

from assay.combine import AssayCombiner
from assay.similarity import SimilarityScorer
from assay.assay2table import AssayToDataFrameConverter
from assay.modseq import stringify_modification   

def score_replicates(assays):
    converter = AssayToDataFrameConverter(columns=[
        {
            'name': 'peptideSequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'modification',
            'path': 'modification',
            'convert': lambda x: str(stringify_modification(x))
        },
        {
            'name': 'glycanStruct',
            'path': 'glycanStruct'
        },
        {
            'name': 'glycanSite',
            'path': 'glycanSite'
        },
        {
            'name': 'precursorCharge',
            'path': 'precursorCharge'
        },
        {
            'name': 'file',
            'path': ['metadata', 'file']
        },
        {
            'name': 'scan',
            'path': ['metadata', 'scan']
        },
        {
            'name': 'rt',
            'path': 'rt'
        }
    ])
    scorer = SimilarityScorer()
    combiner = AssayCombiner()
    
    def score(spectra):
        data = converter.assays_to_dataframe(spectra)
        data.reset_index(drop=True, inplace=True)
        data.sort_values('file', inplace=True)
        spectra = [spectra[i] for i in data.index]
        data.reset_index(drop=True, inplace=True)
    
        data = pd.concat((
            data \
                .iloc[list(itertools.chain.from_iterable((
                    itertools.repeat(i, len(spectra) - i - 1)
                    for i in range(len(spectra) - 1)
                )))] \
                .rename(columns={k: k + '1' for k in ['file', 'scan', 'rt']}) \
                .reset_index(drop=True),
            data \
                .iloc[list(itertools.chain.from_iterable((
                    range(i + 1, len(spectra))
                    for i in range(len(spectra) - 1)
                )))][['file', 'scan', 'rt']] \
                .rename(columns={k: k + '2' for k in ['file', 'scan', 'rt']}) \
                .reset_index(drop=True)
        ), axis=1)
        
        data['delta_rt'] = data['rt2'] - data['rt1']
        
        similarity = scorer.pairwise_similarity(spectra)
        data['intensity_similarity'] = \
            list(itertools.chain.from_iterable(similarity))
        
        return data
    
    return pd.concat((
        score(spectra)
        for spectra in combiner.group_replicates(assays)
    ), axis=0, ignore_index=True)    

# %%
logging.info('scoring replicates')    
    
scores = score_replicates(assays)

logging.info('replicate scoring done') 

# %%
logging.info('saving scoring results: {0}' \
             .format(out_file))

scores.to_csv(out_file, index=False)

logging.info('scoring results saved: {0}' \
             .format(out_file))

logging.info('intensity similarity: median={0}, quantile=({1}, {2})'.format(
    scores['intensity_similarity'].median(),
    scores['intensity_similarity'].quantile(0.25),
    scores['intensity_similarity'].quantile(0.75)
))

logging.info('RT correlation: \n{0}'.format(
    scores.groupby(['file1', 'file2']) \
        .apply(lambda x: x['rt1'].corr(x['rt2']))
))

logging.info('RT difference: IQR={0}, range(95%)={1}'.format(
    scores['delta_rt'].quantile(0.75) - scores['delta_rt'].quantile(0.25),
    scores['delta_rt'].quantile(0.975) - scores['delta_rt'].quantile(0.025),
))
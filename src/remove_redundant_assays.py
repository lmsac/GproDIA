import argparse

parser = argparse.ArgumentParser(
    description='Remove redundant assays.'
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
    '--action', choices=['consensus', 'best', 'first'], default='consensus',
    help='build consensus assays, keep the best replicates, or keep the first occurrence (default: %(default)s)'
)
parser.add_argument(
    '--score', default='totalScore', 
    help='score used for ranking assays (default: %(default)s)'
)
parser.add_argument(
    '--glycan_key', choices=['struct', 'composition'], default='struct',
    help='group replicates by glycan structure or glycan composition (default: %(default)s)'
)

glycansite_group = parser.add_mutually_exclusive_group(required=False)
glycansite_group.add_argument(
    '--use_glycan_site', 
    dest='use_glycan_site', action='store_true', 
    help='consider glycan site (default: %(default)s)'
)
glycansite_group.add_argument(
    '--ignore_glycan_site', 
    dest='use_glycan_site', action='store_false',
    help='ignore glycan site (default: False)'
)
parser.set_defaults(use_glycan_site=True)

withinrun_group = parser.add_mutually_exclusive_group(required=False)
withinrun_group.add_argument(
    '--within_run', 
    dest='within_run', action='store_true',
    help='remove redundant assays within each run (default: %(default)s)'
)
withinrun_group.add_argument(
    '--across_run', 
    dest='within_run', action='store_false',
    help='remove redundant assays across all runs (default: True)'
)
parser.set_defaults(within_run=False)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
action = args.action
score = args.score
glycan_key = args.glycan_key
use_glycan_site = args.use_glycan_site
within_run = args.within_run

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
if globals().get('action', None) is None:
    action = 'consensus'    
    
logging.info('use action: ' + str(action))

if globals().get('score', None) is None:
    score = 'totalScore'

logging.info('use score: ' + str(score))
    
if globals().get('glycan_key', None) is None:
    glycan_key = 'struct'

logging.info('use glycan_key: ' + str(glycan_key))

if globals().get('use_glycan_site', None) is None:
    use_glycan_site = True
    
logging.info('use_glycan_site: ' + str(use_glycan_site))

# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    if action == 'consensus':
        out_file += '_consensus.assay.pickle'
    else:
        out_file += '_nonredundant.assay.pickle'

# %%
from util import save_pickle, load_pickle
from assay.combine import glycopeptide_group_key

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
group_key = glycopeptide_group_key(
    use_glycan_struct=(glycan_key == 'struct'),
    use_glycan_site=use_glycan_site,
    within_run=within_run
)

# %%
if action == 'consensus':
    from assay.consensus import ConsensusAssayCombiner
    combiner = ConsensusAssayCombiner(
        group_key=group_key, 
        replicate_weight=score
    )
elif action == 'best':
    from assay.combine import BestReplicateAssayCombiner
    combiner = BestReplicateAssayCombiner(
        group_key=group_key,
        score=score
    )
else:
    from assay.combine import AssayCombiner
    combiner = AssayCombiner(
        group_key=group_key
    )
    
# %%
if action == 'consensus':
    logging.info('building consensus assays')
else:
    logging.info('removing redundant assays')

assays = combiner.remove_redundant(assays)

logging.info('redundant assays removed: {0} spectra remaining' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)
    
logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

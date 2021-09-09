import argparse

parser = argparse.ArgumentParser(
    description='Filter assays.'
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
    '--swath_windows',
    help='SWATH isolation window file'
)

assay_filter_group = parser.add_argument_group('assay filters') 
assay_filter_group.add_argument(
    '--min_precursor_mz', type=float,
    help='lower m/z limit of precursor ions'
)     
assay_filter_group.add_argument(
    '--max_precursor_mz', type=float,
    help='upper m/z limit of precursor ions'
) 
assay_filter_group.add_argument(
    '--min_fragment_number', default=6, type=int,
    help='remove assays with < N fragments (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--min_peptide_fragment_number', default=3, type=int,
    help='remove assays with < N peptide fragments (default: %(default)s)'
)
assay_filter_group.add_argument(
    '--min_glycan_fragment_number', default=3, type=int,
    help='remove assays with < N glycan fragments (default: %(default)s)'
)

def add_fragment_filter_args(parser):
    for quantify_group in ['main', 'quantify']:
        if quantify_group == 'quantify':
            arg_quantify_group = 'quantify_'
            help_quantify_group = 'quantifying '
        else:
            arg_quantify_group = ''
            help_quantify_group = ''
            
        for prior_group in ['main', 'prior_peptide', 'prior_glycan']:
            if prior_group == 'prior_peptide':
                arg_group = arg_quantify_group + 'prior_peptide_'
                help_group = help_quantify_group + 'prior peptide '
            elif prior_group == 'prior_glycan':
                arg_group = arg_quantify_group + 'prior_glycan_'
                help_group = help_quantify_group + 'prior glycan '
            else:
                arg_group = arg_quantify_group
                help_group = help_quantify_group
            
            frag_filter_group = parser.add_argument_group(help_group + 'fragment filters') 
            if prior_group == 'prior_peptide' or prior_group == 'prior_glycan':
                if quantify_group == 'quantify':
                    default = 6
                else:
                    default = 10
                frag_filter_group.add_argument(
                    '--%sfragment_number' % arg_group, default=default, type=int,
                    help='try to select top N %sfragments' % help_group + ' (default: %(default)s)'
                )
            else:
                if quantify_group == 'quantify':
                    default = 12
                else:
                    default = 20
                frag_filter_group.add_argument(
                    '--%smax_fragment_number' % arg_group, default=default, type=int,
                    help='maximal number of fragments (default: %(default)s)'
                )
            frag_filter_group.add_argument(
                '--%sfragment_type' % arg_group, type=str, nargs='+',
                help='list of %sfragment types' % help_group
            )
            if prior_group == 'main' or prior_group == 'prior_peptide':
                frag_filter_group.add_argument(
                    '--%smin_fragment_amino_acid_number' % arg_group, type=int,
                    help='lower limit of amino acid number of %sfragment ions' % help_group
                )
            frag_filter_group.add_argument(
                '--%sfragment_charge' % arg_group, type=int, nargs='+',
                help='list of allowed charge states of %sfragment ions' % help_group
            )
            frag_filter_group.add_argument(
                '--%sfragment_loss_type' % arg_group, type=str, nargs='+',
                help='list of neutral loss types of %sfragment ions' % help_group
            )  
            frag_filter_group.add_argument(
                '--%smin_fragment_mz' % arg_group, type=float,
                help='lower m/z limit of %sfragment ions' % help_group
            )
            frag_filter_group.add_argument(
                '--%smax_fragment_mz' % arg_group, type=float,
                help='upper m/z limit of %sfragment ions' % help_group
            )    
            if prior_group == 'main' or prior_group == 'prior_glycan':                
                frag_filter_group.add_argument(
                    '--%smin_fragment_monosaccharide_number' % arg_group, type=eval, default=1,
                    help='lower limit of monosaccharide number of %sfragment ions' % help_group + ' (default: %(default)s)'
                )
            if prior_group == 'main':                
                frag_filter_group.add_argument(
                    '--%smin_relative_fragment_intensity' % arg_group, type=float,
                    help='lower relative intensity limit of %sfragment ions' % help_group 
                )

add_fragment_filter_args(parser)

args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out
swath_window_file = args.swath_windows

filter_args = vars(args)
filter_args.pop('in')
filter_args.pop('out')
filter_args.pop('swath_windows')

def arrange_filter_args(filter_args):
    main_args = {
        'prior_peptide_fragment_criteria': {},
        'prior_glycan_fragment_criteria': {},
        'quantifying_transition_criteria': {
            'prior_peptide_fragment_criteria': {},
            'prior_glycan_fragment_criteria': {}
        }
    }
    for k, v in filter_args.items():        
        if k.startswith('quantify_'):
            target = main_args['quantifying_transition_criteria']
            k = k[len('quantify_'):]
        else:
            target = main_args
        if k.startswith('prior_peptide_') and k != 'prior_peptide_fragment_number':                
            k = k[len('prior_peptide_'):]
            target = target['prior_peptide_fragment_criteria']
        elif k.startswith('prior_glycan_') and k != 'prior_glycan_fragment_number':
            k = k[len('prior_glycan_'):]
            target = target['prior_glycan_fragment_criteria']
        
        target[k] = v  
    
    main_args['min_peptide_fragment_criteria'] = main_args['prior_peptide_fragment_criteria']
    main_args['min_glycan_fragment_criteria'] = main_args['prior_glycan_fragment_criteria']
    return main_args
    
filter_criteria = arrange_filter_args(filter_args) 
 
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
    out_file += '_filtered.assay.pickle'

# %%
from util import save_pickle, load_pickle
from assay import GlycoAssayBuilder
import pandas as pd

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
if swath_window_file is not None:
    logging.info('loading SWATH windows: ' + swath_window_file) 
    
    swath_windows = pd.read_csv(swath_window_file, sep='\t')
    
    logging.info('SWATH windows loaded: {0} windows' \
                 .format(len(swath_windows)))
else:
    swath_windows = None
    
# %%
logging.info(
    'filtering assays using the following parameters: \n' + \
    '\n'.join((
        k + '=' + str(v) 
        for k, v in filter_args.items()
        if v is not None
    ))
)

assay_builder = GlycoAssayBuilder()

assays = assay_builder.filter_assays(
    assays,
    swath_windows=swath_windows,
    **filter_criteria
)

logging.info('assays filtered: {0} spectra remaining' \
    .format(len(assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)
    
logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))


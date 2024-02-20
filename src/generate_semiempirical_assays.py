import argparse

parser = argparse.ArgumentParser(
    description='Generate semi-empirical glycopeptide assays.'
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
    '--action', choices=['interchange', 'exchange', 'from_list', 'cross_validation'],
    default='interchange',
    help='interchange peptides and glycans of input assays, ' + \
    'exchange peptide and glycan assays, ' + \
    'generate semi-empirical assays from a list of glycopeptides, ' + \
    'or run cross validation on input assays (default: %(default)s)'
)

parameter_group = parser.add_argument_group('input files for exchange/from-list generation')
parameter_group.add_argument(
    '--peptide_assay', nargs='+',
    help='input peptide assay files'
)
parameter_group.add_argument(
    '--glycan_assay', nargs='+',
    help='input glycan assay files'
)

parameter_group = parser.add_argument_group('input files for from-list generation')
parameter_group.add_argument(
    '--list', nargs='+',
    help='input glycopeptide list files (CSV)'
)

parameter_group = parser.add_argument_group('semi-empirical assay generation parameters')
parameter_group.add_argument(
    '--max_peptide_neighbor_number', default=3, type=int,
    help='maximum peptide neighbor number in kNN-based MS/MS prediction (default: %(default)s)'
)
parameter_group.add_argument(
    '--max_glycan_neighbor_number', default=3, type=int,
    help='maximum glycan neighbor number in kNN-based RT and MS/MS prediction (default: %(default)s)'
)

parameter_group = parser.add_argument_group('parameters for interchange/exchange/cross-validation')
parameter_group.add_argument(
    '--min_peptide_occurrence', default=3, type=int,
    help='ignore assays with peptide occurrence < N (default: %(default)s)'
)
parameter_group.add_argument(
    '--min_glycan_occurrence', default=3, type=int,
    help='ignore assays with glycan occurrence < N (default: %(default)s)'
)
parameter_group.add_argument(
    '--top_n_assays_by_occurrence', default=10000, type=int,
    help='for interchange/exchange, select top N assays by peptide and glycan occurrence (default: %(default)s)'
)
parameter_group.add_argument(
    '--random_select_n_assays', default=10000, type=int,
    help='for interchange/exchange/cross-validation, randomly select N assays (default: %(default)s)'
)

parser.add_argument(
    '--glycan_key', choices=['struct', 'composition'], default='struct',
    help='consider glycan as glycan structure or glycan composition (default: %(default)s)'
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


args = parser.parse_args()
assay_files = getattr(args, 'in')
out_file = args.out

interchange = args.action == 'interchange'
exchange = args.action == 'exchange'
from_list = args.action == 'from_list'
cross_validation = args.action == 'cross_validation'

peptide_assay_files = args.peptide_assay
glycan_assay_files = args.glycan_assay

glycopeptide_list_files = args.list

max_peptide_neighbor_number = args.max_peptide_neighbor_number
max_glycan_neighbor_number = args.max_glycan_neighbor_number

min_peptide_occurrence = args.min_peptide_occurrence
min_glycan_occurrence = args.min_glycan_occurrence
top_n_assays_by_occurrence = args.top_n_assays_by_occurrence
random_select_n_assays = args.random_select_n_assays

glycan_key = args.glycan_key
use_glycan_site = args.use_glycan_site

# %%
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%
from util import list_files

if interchange or cross_validation:
    if globals().get('assay_files', None) is None:
        assay_files = list_files(
            path='.',
            pattern='\\.assay\\.pickle$'
        )

    if len(assay_files) == 0:
        raise ValueError('no assay files')

if exchange:
    if len(peptide_assay_files) == 0:
        raise ValueError('no peptide assay files')
    if len(glycan_assay_files) == 0:
        raise ValueError('no glycan assay files')

if from_list:
    if len(glycopeptide_list_files) == 0:
        raise ValueError('no glycopeptide list files')

# %%
import os

if interchange or cross_validation:
    if globals().get('out_file', None) is None:
        out_file = os.path.splitext(assay_files[0])[0]
        if out_file.endswith('.assay'):
            out_file = out_file[:-len('.assay')]
        if len(assay_files) > 1:
            out_file += '_' + str(len(assay_files))
        if cross_validation:
            out_file += '_semiempirical_validation.assay.pickle'
        else:
            out_file += '_semiempirical.assay.pickle'

if exchange:
    if globals().get('out_file', None) is None:
        out_file = os.path.splitext(peptide_assay_files[0])[0]
        if out_file.endswith('.assay'):
            out_file = out_file[:-len('.assay')]
        if len(peptide_assay_files) > 1:
            out_file += '_' + str(len(peptide_assay_files))
        out_file += '_' + \
            os.path.splitext(os.path.basename(glycan_assay_files[0]))[0]
        if out_file.endswith('.assay'):
            out_file = out_file[:-len('.assay')]
        if len(glycan_assay_files) > 1:
            out_file += '_' + str(len(glycan_assay_files))
        out_file += '_semiempirical.assay.pickle'

if from_list:
    if globals().get('out_file', None) is None:
        out_file = os.path.splitext(glycopeptide_list_files[0])[0]
        if out_file.endswith('.glycopeptide'):
            out_file = out_file[:-len('.glycopeptide')]
        if len(glycopeptide_list_files) > 1:
            out_file += '_' + str(len(glycopeptide_list_files))
        out_file += '_semiempirical.assay.pickle'

# %%
from util import save_pickle, load_pickle

# %%
if interchange or cross_validation or from_list:
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
if exchange or from_list:
    if peptide_assay_files is not None:
        peptide_assays = []
        for assay_file in peptide_assay_files:
            logging.info('loading peptide assays: ' + assay_file)

            assay_data = load_pickle(assay_file)
            peptide_assays.extend(assay_data)

            logging.info('peptide assays loaded: {0}, {1} spectra' \
                .format(assay_file, len(assay_data)))

        logging.info('peptide assays loaded: {0} spectra totally' \
            .format(len(peptide_assays)))
    else:
        peptide_assays = None

    if glycan_assay_files is not None:
        glycan_assays = []
        for assay_file in glycan_assay_files:
            logging.info('loading glycan assays: ' + assay_file)

            assay_data = load_pickle(assay_file)
            glycan_assays.extend(assay_data)

            logging.info('glycan assays loaded: {0}, {1} spectra' \
                .format(assay_file, len(assay_data)))

        logging.info('glycan assays loaded: {0} spectra totally' \
            .format(len(glycan_assays)))
    else:
        glycan_assays = None

# %%
if from_list:
    import pandas as pd

    logging.info('loading glycopeptide list file(s): ' + \
                 '; '.join(glycopeptide_list_files))

    glycopeptide_data = pd.concat(
        (pd.read_csv(f) for f in glycopeptide_list_files),
        ignore_index=True
    )

    logging.info('glycopeptide list file(s) loaded: {0} glycopeptides' \
        .format(len(glycopeptide_data)))


# %%
if interchange:
    from assay.semiempirical import interchange_peptide_glycan

    logging.info('generating glycopeptides by interchanging peptides and glycans')

    new_assays, glycopeptide_table  = interchange_peptide_glycan(
        assays=assays,
        return_generator=True,
        return_glycopeptide_table=True,
        min_peptide_occurrence=min_peptide_occurrence,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=(glycan_key == 'struct'),
        use_glycan_site=use_glycan_site,
        max_peptide_neighbor_number=max_peptide_neighbor_number,
        max_glycan_neighbor_number=max_glycan_neighbor_number,
        top_n_assays_by_occurrence=top_n_assays_by_occurrence,
        random_select_n_assays=random_select_n_assays
    )

    logging.info('generating semi-empirical assays of {0} glycopeptides' \
                 .format(len(glycopeptide_table)))

    new_assays = list(new_assays)

    logging.info('semi-empirical assays generated: {0} spectra' \
        .format(len(new_assays)))

# %%
if exchange:
    from assay.semiempirical import exchange_peptide_glycan

    logging.info('generating by exchanging peptides and glycans')

    new_assays, glycopeptide_table = exchange_peptide_glycan(
        peptide_assays=peptide_assays,
        glycan_assays=glycan_assays,
        return_generator=True,
        return_glycopeptide_table=True,
        min_peptide_occurrence=min_peptide_occurrence,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=(glycan_key == 'struct'),
        use_glycan_site=use_glycan_site,
        max_peptide_neighbor_number=max_peptide_neighbor_number,
        max_glycan_neighbor_number=max_glycan_neighbor_number,
        top_n_assays_by_occurrence=top_n_assays_by_occurrence,
        random_select_n_assays=random_select_n_assays
    )

    logging.info('generating semi-empirical assays of {0} glycopeptides' \
                 .format(len(glycopeptide_table)))

    import tqdm
    new_assays = list(tqdm.tqdm(new_assays, total=len(glycopeptide_table)))

    logging.info('semi-empirical assays generated: {0} spectra' \
        .format(len(new_assays)))

# %%
if from_list:
    from assay.semiempirical import generate_semiempirical_assays

    logging.info('generating semi-empirical assays from a list of glycopeptides')

    new_assays, glycopeptide_table = generate_semiempirical_assays(
        data=glycopeptide_data,
        assays=assays,
        peptide_assays=peptide_assays,
        glycan_assays=glycan_assays,
        return_generator=True,
        return_glycopeptide_table=True,
        min_peptide_occurrence=min_peptide_occurrence,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=(glycan_key == 'struct'),
        use_glycan_site=use_glycan_site,
        max_peptide_neighbor_number=max_peptide_neighbor_number,
        max_glycan_neighbor_number=max_glycan_neighbor_number
    )

    logging.info('generating semi-empirical assays of {0} glycopeptides' \
                 .format(len(glycopeptide_table)))

    new_assays = list(new_assays)

    logging.info('semi-empirical assays generated: {0} spectra' \
        .format(len(new_assays)))

# %%
if cross_validation:
    from assay.semiempirical import \
        generate_semiempirical_assays_for_cross_validation

    logging.info('generating semi-empirical assays for cross validation')

    new_assays, glycopeptide_table = generate_semiempirical_assays_for_cross_validation(
        assays=assays,
        include_index=True,
        return_generator=True,
        return_glycopeptide_table=True,
        min_peptide_occurrence=min_peptide_occurrence,
        min_glycan_occurrence=min_glycan_occurrence,
        max_peptide_neighbor_number=max_peptide_neighbor_number,
        max_glycan_neighbor_number=max_glycan_neighbor_number,
        random_select_n_assays=random_select_n_assays
    )

    logging.info('generating semi-empirical assays of {0} glycopeptides' \
                 .format(len(glycopeptide_table)))

    new_assays = list(new_assays)

    logging.info('semi-empirical assays generated: {0} spectra' \
        .format(len(new_assays)))

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

if not cross_validation:
    save_pickle(new_assays, out_file)
else:
    save_pickle([t[1] for t in new_assays], out_file)

logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(new_assays)))

# %%
if cross_validation:
    from assay.similarity import SimilarityScorer
    import pandas as pd

    logging.info('running cross validation')

    scorer = SimilarityScorer()
    scores = pd.DataFrame.from_records(
        [
            (
                i, t[0],
                scorer.similarity(t[1], assays[t[0]]),
                t[1]['rt'], assays[t[0]]['rt'], t[1]['rt'] - assays[t[0]]['rt']
            )
            for i, t in enumerate(new_assays)
        ],
        columns=[
            'index_semiempirical', 'index_empirical',
            'intensity_similarity',
            'rt_semiempirical', 'rt_empirical', 'delta_rt'
        ]
    )

    logging.info('cross validation done')

    scores = pd.merge(
        glycopeptide_table,
        scores,
        left_on='index', right_on='index_empirical',
        how='right'
    ).drop(columns=['index'])

    out_score_file = os.path.splitext(out_file)[0]
    if out_score_file.endswith('.assay'):
        out_score_file = out_score_file[:-len('.assay')]
    out_score_file += '.score.csv'

    logging.info('saving cross validation results: {0}' \
                 .format(out_score_file))

    scores.to_csv(out_score_file, index=False)

    logging.info('cross validation results saved: {0}' \
                 .format(out_score_file))

    logging.info('intensity similarity: median={0}, quantile=({1}, {2})'.format(
        scores['intensity_similarity'].median(),
        scores['intensity_similarity'].quantile(0.25),
        scores['intensity_similarity'].quantile(0.75)
    ))

    logging.info('RT correlation: pearson={0}'.format(
        scores['rt_semiempirical'].corr(scores['rt_empirical'])
    ))

    logging.info('RT difference: IQR={0}, range(95%)={1}'.format(
        scores['delta_rt'].quantile(0.75) - scores['delta_rt'].quantile(0.25),
        scores['delta_rt'].quantile(0.975) - scores['delta_rt'].quantile(0.025),
    ))

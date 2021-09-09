import pandas as pd
import numpy as np
import sqlite3
import os
import click
from scipy.interpolate import make_interp_spline

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

from spectra.mzmlreader import MzmlReader
from assay.annotation import match_fragments, match_peaks


def read_feature_transition(infile, run_id,
                            feature_id=None,
                            transition_group_id=None,
                            glycoform=False,
                            max_peakgroup_rank=None,
                            max_rs_peakgroup_qvalue=None,
                            max_glycoform_qvalue=None,
                            include_decoy=False):
    con = sqlite3.connect(infile)

    query_column = '''
SELECT RUN.ID AS run_id,
       RUN.FILENAME AS filename,
       PRECURSOR.ID AS transition_group_id,

       PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
       PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
       GLYCAN.GLYCAN_STRUCT AS GlycanStruct,
       GLYCAN.GLYCAN_COMPOSITION AS GlycanComposition,
       GLYCOPEPTIDE.GLYCAN_SITE AS GlycanSite,
       PRECURSOR.CHARGE AS Charge,
       PRECURSOR.PRECURSOR_MZ AS mz,
       PRECURSOR.LIBRARY_RT AS LibraryRT,

       FEATURE.ID AS id,
       FEATURE.EXP_RT AS RT,
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       SCORE_MS2.RANK AS peak_group_rank,

       TRANSITION.ID AS transition_id,
       TRANSITION.ANNOTATION AS FragmentAnnotation,
       TRANSITION.PRODUCT_MZ AS ProductMZ,
       TRANSITION.LIBRARY_INTENSITY AS LibraryIntensity,
       TRANSITION.QUANTIFYING AS quantifyingTransition
'''

    query_table = '''
FROM PRECURSOR
INNER JOIN PRECURSOR_GLYCOPEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_GLYCOPEPTIDE_MAPPING.PRECURSOR_ID
INNER JOIN GLYCOPEPTIDE ON PRECURSOR_GLYCOPEPTIDE_MAPPING.GLYCOPEPTIDE_ID = GLYCOPEPTIDE.ID
INNER JOIN GLYCOPEPTIDE_PEPTIDE_MAPPING ON GLYCOPEPTIDE.ID = GLYCOPEPTIDE_PEPTIDE_MAPPING.GLYCOPEPTIDE_ID
INNER JOIN PEPTIDE ON GLYCOPEPTIDE_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
INNER JOIN GLYCOPEPTIDE_GLYCAN_MAPPING ON GLYCOPEPTIDE.ID = GLYCOPEPTIDE_GLYCAN_MAPPING.GLYCOPEPTIDE_ID
INNER JOIN GLYCAN ON GLYCOPEPTIDE_GLYCAN_MAPPING.GLYCAN_ID = GLYCAN.ID

INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN TRANSITION ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID
INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
'''

    conditions = {
        'RUN.ID': run_id,
        'FEATURE.ID': feature_id,
        'PRECURSOR.ID': transition_group_id,
        'SCORE_MS2.RANK <=': max_peakgroup_rank,
        'SCORE_MS2.QVALUE <': max_rs_peakgroup_qvalue
    }
    if not include_decoy:
        conditions['PRECURSOR.DECOY'] = 0
        conditions['TRANSITION.DECOY'] = 0

    if glycoform:
        query_column += ''',
           TRANSITION.IDENTIFYING AS identifyingTransition,
           IFNULL(GLYCOFORM.GLYCOPEPTIDE_ID, -1) AS glycoform,
           SCORE_MS2.QVALUE AS ms2_m_score,
           SCORE_GLYCOFORM.QVALUE AS m_score
'''

        query_table += '''
INNER JOIN SCORE_GLYCOFORM ON SCORE_GLYCOFORM.FEATURE_ID = FEATURE.ID
                          AND SCORE_GLYCOFORM.GLYCOPEPTIDE_ID = GLYCOPEPTIDE.ID
LEFT JOIN FEATURE_TRANSITION ON TRANSITION.ID = FEATURE_TRANSITION.TRANSITION_ID
                            AND FEATURE.ID = FEATURE_TRANSITION.FEATURE_ID
LEFT JOIN (
    SELECT TRANSITION_ID,
           GLYCOPEPTIDE_ID
    FROM TRANSITION_GLYCOPEPTIDE_MAPPING
    INNER JOIN GLYCOPEPTIDE AS GLYCOPEPTIDE ON GLYCOPEPTIDE.ID = TRANSITION_GLYCOPEPTIDE_MAPPING.GLYCOPEPTIDE_ID
) AS GLYCOFORM ON GLYCOFORM.TRANSITION_ID = TRANSITION.ID
              AND GLYCOFORM.GLYCOPEPTIDE_ID = GLYCOPEPTIDE.ID
'''

        # conditions['TRANSITION.IDENTIFYING'] = 1
        conditions['SCORE_GLYCOFORM.QVALUE <'] = max_glycoform_qvalue

        conditions['TRANSITION.DETECTING = 1 OR IFNULL(FEATURE_TRANSITION.TRANSITION_ID, -1) !='] = -1
    else:
        query_column += ''',
           SCORE_MS2.QVALUE AS m_score
'''

        conditions['TRANSITION.DETECTING'] = 1

    where = [
        k + ' ' + str(v) if \
            k.endswith('=') or k.endswith('<') or k.endswith('>') \
        else \
        k + ' IN (' + ','.join(map(str, v)) + ')' if \
            isinstance(v, list) or \
            isinstance(v, tuple) or \
            isinstance(v, set) \
        else \
        k + ' = ' + str(v)
        for k, v in conditions.items()
        if v is not None
    ]
    query_where = 'WHERE (' + ')\n AND ('.join(where) + ')\n'

    query_order = '''
ORDER BY transition_group_id,
         run_id,
         peak_group_rank;
'''

    query = query_column + query_table + query_where + query_order
    data = pd.read_sql_query(query, con)
    return data


def extract_feature_chromatograms(spectra, feature_data, swath_windows,
                                  tolerance=20, tolerance_unit='ppm',
                                  ms1_tolerance=10, ms1_tolerance_unit='ppm',
                                  criteria='mostintense'):
    def to_assay(data):
        result = {
            'peptideSequence': data['Sequence'].values[0],
            'modifiedSequence': data['FullPeptideName'].values[0],
            'glycanStruct': data['GlycanStruct'].values[0],
            'glycanSite': int(data['GlycanSite'].values[0]),
            'precursorCharge': int(data['Charge'].values[0]),
            'precursorMZ': float(data['mz'].values[0]),
            'rt': float(data['RT'].values[0]),

            'metadata': {
                'file': data['filename'].values[0],
                'runID': int(data['run_id'].values[0]),
                'transitionGroupID': int(data['transition_group_id'].values[0]),
                'qvalue': float(data['m_score'].values[0]),
                'leftWidth': float(data['leftWidth'].values[0]),
                'rightWidth': float(data['rightWidth'].values[0]),
                'libraryRT': float(data['LibraryRT'].values[0]),
                'libraryIntensity': data['LibraryIntensity'] \
                    .values.tolist(),
            },
            'fragments' : {
                'fragmentMZ': data['ProductMZ'].values.tolist(),
                'fragmentAnnotation': data['FragmentAnnotation'] \
                    .values.tolist(),
                'quantifyingTransition': data['quantifyingTransition'] \
                    .values.tolist()
            }
        }
        if 'identifyingTransition' in data.columns:
            result['fragments']['identifyingTransition'] = \
                data['identifyingTransition'].values.tolist()
        if 'glycoform' in data.columns:
            result['fragments']['glycoform'] = \
                data['glycoform'].values.tolist()
        if 'ms2_m_score' in data.columns:
            result['metadata']['ms2_qvalue'] = \
                float(data['ms2_m_score'].values[0])

        swath_window_index = \
            np.max(np.append(np.where(
                (swath_windows['start'] < result['precursorMZ']) & \
                (swath_windows['end'] > result['precursorMZ'])
            )[0], -1))

        result['metadata'].update({
            'isolationWindowLowerLimit': \
                float(swath_windows['start'][swath_window_index]),
            'isolationWindowUpperLimit': \
                float(swath_windows['end'][swath_window_index])
        })
        return result


    def extract_ms2(features, spec, rt, rt_margin=0.2, min_rt_window=60):
        for feature in features:
            rt_left = feature['metadata']['leftWidth']
            rt_right = feature['metadata']['rightWidth']

            if rt_left > feature['rt']:
                rt_left = feature['rt']
            if rt_right < feature['rt']:
                rt_right = feature['rt']

            if rt_margin is not None:
                rt_left -= (feature['metadata']['rightWidth'] - \
                     feature['metadata']['leftWidth']) * rt_margin
                rt_right += (feature['metadata']['rightWidth'] - \
                     feature['metadata']['leftWidth']) * rt_margin

            if min_rt_window is not None and rt_right - rt_left < min_rt_window:
                rt_left = feature['rt'] - (feature['rt'] - rt_left) * \
                    min_rt_window / (rt_right - rt_left)
                rt_right = feature['rt'] + (rt_right - feature['rt']) * \
                    min_rt_window / (rt_right - rt_left)

            if rt_left > rt or rt_right < rt:
                continue

            if feature['metadata']['isolationWindowLowerLimit'] \
                > spec['metadata']['isolationWindowTargetMZ'] or \
                feature['metadata']['isolationWindowUpperLimit'] \
                < spec['metadata']['isolationWindowTargetMZ']:
                    continue

            chromatograms = feature.get('chromatograms', None)
            if chromatograms is None:
                chromatograms = {
                    'rt': [],
                    'intensity': list((
                        [] for x in feature['fragments']['fragmentMZ']
                    )),
                    'mz' : list((
                        [] for x in feature['fragments']['fragmentMZ']
                    ))
                }
                feature['chromatograms'] = chromatograms

            chromatograms['rt'].append(rt)
            for x in chromatograms['intensity']:
                x.append(0)
            for x in chromatograms['mz']:
                x.append(None)

            index = match_fragments(
                feature, spec,
                tolerance=tolerance, tolerance_unit=tolerance_unit,
                criteria=criteria,
                all1=True, all2=False
            )
            for i, j in index:
                if j is not None:
                    chromatograms['intensity'][i][-1] = \
                        spec['fragments']['fragmentIntensity'][j]
                    chromatograms['mz'][i][-1] = \
                        spec['fragments']['fragmentMZ'][j]


    def extract_ms1(features, spec, rt):
        for feature in features:
            if feature['metadata']['leftWidth'] > rt or \
                feature['metadata']['rightWidth'] < rt:
                continue

            chromatograms = feature.get('ms1Chromatogram', None)
            if chromatograms is None:
                chromatograms = {
                    'rt': [],
                    'intensity': [[]],
                    'mz' : [[]]
                }
                feature['ms1Chromatogram'] = chromatograms

            chromatograms['rt'].append(rt)
            for x in chromatograms['intensity']:
                x.append(0)
            for x in chromatograms['mz']:
                x.append(None)

            index = match_peaks(
                { 'peaks': { 'mz': [feature['precursorMZ']] } },
                spec,
                tolerance=ms1_tolerance, tolerance_unit=ms1_tolerance_unit,
                criteria=criteria,
                all1=True, all2=False
            )
            for i, j in index:
                if j is not None:
                    chromatograms['intensity'][i][-1] = \
                        spec['peaks']['intensity'][j]
                    chromatograms['mz'][i][-1] = \
                        spec['peaks']['mz'][j]


    features = [
        to_assay(data)
        for id, data in feature_data.groupby(by=['id'], sort=False)
    ]

    rt_start = min((
        feature['metadata']['leftWidth']
        for feature in features
    ))
    rt_end = max((
        feature['metadata']['rightWidth']
        for feature in features
    ))

    for spec in spectra:
        rt = spec['rt'] * 60
        if rt < rt_start:
            continue
        if rt > rt_end:
            break

        if spec['msLevel'] == 2:
            extract_ms2(features, spec, rt)
        elif spec['msLevel'] == 1:
            extract_ms1(features, spec, rt)
        else:
            continue

    for feature in features:
        feature['fragments']['fragmentIntensity'] = \
            [sum(x) for x in feature['chromatograms']['intensity']]
        feature['metadata']['libraryProductMZ'] = \
            feature['fragments']['fragmentMZ']
        feature['fragments']['fragmentMZ'] = \
            [
                float(
                    np.mean([y for y in x if y is not None]) \
                    if any(y is not None for y in x) \
                    else np.nan
                )
                for x in feature['chromatograms']['mz']
            ]

        feature['precursorIntensity'] = \
            sum(feature['ms1Chromatogram']['intensity'][0])
        feature['metadata']['libraryPrecursorMZ'] = \
            feature['precursorMZ']
        feature['precursorMZ'] = float(
            np.mean([
                y for y in feature['ms1Chromatogram']['mz'][0]
                if y is not None
            ]) \
            if any(y is not None for y in feature['ms1Chromatogram']['mz'][0]) \
            else np.nan
        )

    return features


def plot_feature_ms2_chromatograms(ax, feature,
                                   identifying_transitions=True,
                                   quantifying_transitions=True,
                                   nonquantifying_detecting_transitions=True,
                                   other_glycoforms_transitions=False,
                                   including_ms1=False,
                                   title=None, xlabel=None, ylabel=None,
                                   legend=True, smooth=True,
                                   **kwargs):
    if title is None:
        title = str(feature['modifiedSequence']) + ' / ' + \
            str(feature['precursorCharge']) + '\n' + \
            str(feature['glycanStruct']) + ' @ ' + \
            str(feature['glycanSite']) + '\n'
    if xlabel is None:
        xlabel = 'Retention Time'
    if ylabel is None:
        ylabel = 'Intensity'

    if ax is None:
        fig, ax = plt.subplots()

    rt = feature['chromatograms']['rt']
    if smooth and len(rt) >= 3:
        rt_smooth = np.linspace(np.min(rt), np.max(rt), 50)

    quantifying = feature['fragments'].get('quantifyingTransition', None)
    identifying = feature['fragments'].get('identifyingTransition', None)
    glycoform = feature['fragments'].get('glycoform', None)

    for i, annot in enumerate(feature['fragments']['fragmentAnnotation']):
        if feature['fragments']['fragmentIntensity'][i] == 0:
           continue

        if identifying is not None and identifying[i] == 1:
            if not identifying_transitions:
                continue
            if glycoform is not None and glycoform[i] == -1:
                if not other_glycoforms_transitions:
                    continue
                linestyle = '--'
            else:
                linestyle = '-'
            linewidth = 0.75

        elif quantifying is not None and quantifying[i] == 1:
            if not quantifying_transitions:
                continue
            linestyle = '-'
            linewidth = 1.5

        else:
            if not nonquantifying_detecting_transitions:
                continue
            linestyle = '-'
            linewidth = 0.75


        intensity = feature['chromatograms']['intensity'][i]

        if smooth and len(rt) >= 3:
            intensity_smooth = make_interp_spline(rt, intensity, 2)(rt_smooth)
            intensity_smooth = np.maximum(intensity_smooth, 0.0)
            ax.plot(
                rt_smooth, intensity_smooth, label=annot,
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )
        else:
            ax.plot(
                rt, intensity, label=annot,
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )

    if including_ms1:
        rt = feature['ms1Chromatogram']['rt']
        if smooth and len(rt) >= 3:
            rt_smooth = np.linspace(np.min(rt), np.max(rt), 50)

        intensity = feature['ms1Chromatogram']['intensity'][0]

        color = 'red'
        annot = 'MS1'

        if smooth and len(rt) >= 3:
            intensity_smooth = make_interp_spline(rt, intensity, 2)(rt_smooth)
            intensity_smooth = np.maximum(intensity_smooth, 0.0)
            ax.plot(
                rt_smooth, intensity_smooth, label=annot, color=color,
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )
        else:
            ax.plot(
                rt, intensity, label=annot, color=color,
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )


    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    if isinstance(legend, dict):
        ax.legend(**legend)
    elif legend:
        ax.legend()

    return ax


def plot_feature_ms2_spectrum(ax, feature,
                              identifying_transitions=True,
                              quantifying_transitions=True,
                              nonquantifying_detecting_transitions=True,
                              other_glycoforms_transitions=False,
                              title=None, xlabel=None, ylabel=None,
                              **kwargs):
    def normalize(x):
        m = max([y if y is not None else 0 for y in x])
        if m > 0:
            x = [y / m if y is not None else None for y in x]
        return x

    if title is None:
        title = str(feature['modifiedSequence']) + '/' + \
            str(feature['precursorCharge']) + '\n' + \
            str(feature['glycanSite']) + ',' + \
            str(feature['glycanStruct']) + '\n'
    if xlabel is None:
        xlabel = 'm/z'
    if ylabel is None:
        ylabel = 'Relative Intensity'

    if ax is None:
        fig, ax = plt.subplots()

    mz = feature['fragments']['fragmentMZ']
    intensity = normalize(feature['fragments']['fragmentIntensity'])

    quantifying = feature['fragments'].get('quantifyingTransition', None)
    identifying = feature['fragments'].get('identifyingTransition', None)
    glycoform = feature['fragments'].get('glycoform', None)

    for i, annot in enumerate(feature['fragments']['fragmentAnnotation']):
        if not (mz[i] > 0 and intensity[i] > 0):
            continue

        if identifying is not None and identifying[i] == 1:
            if not identifying_transitions:
                continue
            if glycoform is not None and glycoform[i] == -1:
                if not other_glycoforms_transitions:
                    continue
                linestyle = '--'
                color = 'r'
            else:
                linestyle = '-'
                color = 'g'
            linewidth = 0.75

        elif quantifying is not None and quantifying[i] == 1:
            if not quantifying_transitions:
                continue
            color = 'b'
            linestyle = '-'
            linewidth = 1.5

        else:
            if not nonquantifying_detecting_transitions:
                continue
            color = 'b'
            linestyle = '-'
            linewidth = 0.75

        ax.vlines(
            x=mz[i],
            ymin=0, ymax=intensity[i],
            color=color,
            linewidth=linewidth, linestyle=linestyle,
            **kwargs
        )

        ax.text(
            mz[i], intensity[i],
            annot, fontsize=5,
            ha ='center', va ='bottom',
            rotation=90
        )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yticks(
        np.linspace(0, 1, 6),
        list(map(str, list(range(100, 0, -25)) + list(range(0, 101, 25))))
    )

    return ax


def plot_feature_ms2_spectrum_comparison(ax, feature,
                                         identifying_transitions=True,
                                         quantifying_transitions=True,
                                         nonquantifying_detecting_transitions=True,
                                         other_glycoforms_transitions=False,
                                         title=None, xlabel=None, ylabel=None,
                                         **kwargs):
    def normalize(x):
        m = max([y if y is not None else 0 for y in x])
        if m > 0:
            x = [y / m if y is not None else None for y in x]
        return x

    if title is None:
        title = str(feature['modifiedSequence']) + '/' + \
            str(feature['precursorCharge']) + '\n' + \
            str(feature['glycanSite']) + ',' + \
            str(feature['glycanStruct']) + '\n'
    if xlabel is None:
        xlabel = 'm/z'
    if ylabel is None:
        ylabel = 'Relative Intensity'

    if ax is None:
        fig, ax = plt.subplots()

    mz = feature['fragments']['fragmentMZ']
    intensity = normalize(feature['fragments']['fragmentIntensity'])
    lib_mz = feature['metadata']['libraryProductMZ']
    lib_intensity = normalize(feature['metadata']['libraryIntensity'])

    quantifying = feature['fragments'].get('quantifyingTransition', None)
    identifying = feature['fragments'].get('identifyingTransition', None)
    glycoform = feature['fragments'].get('glycoform', None)

    for i, annot in enumerate(feature['fragments']['fragmentAnnotation']):
        if identifying is not None and identifying[i] == 1:
            if not identifying_transitions:
                continue
            if glycoform is not None and glycoform[i] == -1:
                if not other_glycoforms_transitions:
                    continue
                linestyle = '--'
                color = 'y'
            else:
                linestyle = '-'
                color = 'g'
            linewidth = 0.75

        elif quantifying is not None and quantifying[i] == 1:
            if not quantifying_transitions:
                continue
            color = 'b'
            linestyle = '-'
            linewidth = 1.5

        else:
            if not nonquantifying_detecting_transitions:
                continue
            color = 'b'
            linestyle = '-'
            linewidth = 0.75

        if mz[i] > 0 and intensity[i] > 0:
            ax.vlines(
                x=mz[i],
                ymin=0, ymax=intensity[i],
                color=color,
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )

            if lib_intensity[i] is None or \
                not lib_intensity[i] > 0 or \
                intensity[i] > lib_intensity[i]:
                ax.text(
                    mz[i], intensity[i],
                    annot, fontsize=5,
                    ha ='center', va ='bottom',
                    rotation=90
                )

        if lib_mz[i] > 0 and \
            lib_intensity[i] is not None and \
            lib_intensity[i] > 0:
            ax.vlines(
                x=lib_mz[i],
                ymin=0, ymax=-lib_intensity[i],
                color='r',
                linewidth=linewidth, linestyle=linestyle,
                **kwargs
            )

            if intensity[i] is None or \
                not intensity[i] > 0 or \
                intensity[i] <= lib_intensity[i]:
                ax.text(
                    lib_mz[i], -lib_intensity[i],
                    annot, fontsize=5,
                    ha ='center', va ='top',
                    rotation=90
                )

    ax.axhline(y=0, color='grey')

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yticks(
        np.linspace(-1, 1, 9),
        list(map(str, list(range(100, 0, -25)) + list(range(0, 101, 25))))
    )

    return ax


def plot_peakgroup_features(features, pdf_path, glycoform=False):
    if plt is None:
        raise ImportError("Error: The matplotlib package is required to create a report.")

    with PdfPages(pdf_path) as pdf:
        for idx, feature in enumerate(features):
            if idx % 3 == 0:
                fig = plt.figure(figsize=(10, 15))
                fig.subplots_adjust(hspace=.75)

            ax = fig.add_subplot(3, 2, 1 + idx % 3 * 2)
            if not glycoform:
                plot_feature_ms2_chromatograms(
                    ax, feature,
                    legend=dict(fontsize=5),
                    quantifying_transitions=True,
                    identifying_transitions=False,
                    nonquantifying_detecting_transitions=True
                )
            else:
                plot_feature_ms2_chromatograms(
                    ax, feature,
                    legend=dict(fontsize=5),
                    quantifying_transitions=1,
                    identifying_transitions=True,
                    nonquantifying_detecting_transitions=False
                )

            ax = fig.add_subplot(3, 2, 1 + idx % 3 * 2 + 1)
            if not glycoform:
                plot_feature_ms2_spectrum_comparison(
                    ax, feature,
                    title= \
                        'TransitionGroup: ' + str(feature['metadata']['transitionGroupID']) + '\n' + \
                        'Q-value: ' + str(feature['metadata']['qvalue']) + '\n',
                    quantifying_transitions=True,
                    identifying_transitions=False,
                    nonquantifying_detecting_transitions=True
                )
            else:
                plot_feature_ms2_spectrum(
                    ax, feature,
                    title= \
                        'TransitionGroup: ' + str(feature['metadata']['transitionGroupID']) + '\n' + \
                        'Q-value: ' + str(feature['metadata']['qvalue']) + '\n',
                    quantifying_transitions=False,
                    identifying_transitions=True,
                    nonquantifying_detecting_transitions=False
                )

            if idx % 3 == 2 or idx == len(features) - 1:
                pdf.savefig(fig)
                plt.close(fig)


def export_peakgroup_plots(osw_file, mzml_file, swath_window_file,
                           transition_group_file,
                           glycoform=False,
                           max_peakgroup_rank=None,
                           max_rs_peakgroup_qvalue=None,
                           max_glycoform_qvalue=None,
                           include_decoy=False,
                           tolerance=20, tolerance_unit='ppm'):
    swath_windows = pd.read_csv(swath_window_file, sep='\t')

    run_id = None
    feature_id = None
    transition_group_id = None

    if transition_group_file is not None:
        result = pd.read_csv(transition_group_file, sep='\t')
        if 'filename' in result.columns:
            result = result.loc[
                result['filename'].str.endswith(os.path.basename(mzml_file))
            ]
            if 'run_id' in result.columns and result['run_id'].dtype.kind == 'i':
                run_id = result['run_id'].iloc[0]

        if 'decoy' in result.columns and not include_decoy:
            result = result.loc[
                result['decoy'] == 0
            ]
        if 'peak_group_rank' in result.columns and max_peakgroup_rank is not None:
            result = result.loc[
                result['peak_group_rank'] <= max_peakgroup_rank
            ]
        if 'id' in result.columns:
            feature_id = list(result['id'].drop_duplicates())
        else:
            transition_group_id = \
                list(result['transition_group_id'].drop_duplicates())

    feature_data = read_feature_transition(
        osw_file,
        run_id=run_id,
        feature_id=feature_id,
        transition_group_id=transition_group_id,
        glycoform=glycoform,
        max_peakgroup_rank=max_peakgroup_rank,
        max_rs_peakgroup_qvalue=max_rs_peakgroup_qvalue,
        max_glycoform_qvalue=max_glycoform_qvalue,
        include_decoy=include_decoy
    )
    if run_id is None:
        feature_data = feature_data.loc[
            feature_data['filename'].str.endswith(os.path.basename(mzml_file))
        ]

    def mzml_as_iterable(path):
        reader = MzmlReader(path)
        i = 0
        click.echo("Info: Scan", nl=False)
        while True:
            spec = reader.read_spectrum()
            if spec is None:
                break

            i += 1
            if i % 1000 == 0:
                click.echo(" {0}".format(i), nl=False)
            yield spec

        click.echo("")
        click.echo("Info: {0} loaded, {1} scans".format(path, i))

    spectra = mzml_as_iterable(mzml_file)

    features = extract_feature_chromatograms(
        spectra, feature_data, swath_windows,
        tolerance=tolerance, tolerance_unit=tolerance_unit
    )

    outfile = os.path.basename(mzml_file).split('.mzML')[0]
    if transition_group_file is not None:
        outfile += '_' + os.path.basename(transition_group_file).split('.tsv')[0]
    outfile += '_peakgroup_plots.pdf'
    plot_peakgroup_features(features, outfile, glycoform=glycoform)
    click.echo("Info: {0} saved.".format(outfile))

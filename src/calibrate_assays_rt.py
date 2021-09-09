import argparse

parser = argparse.ArgumentParser(
    description='Calibrate retention time.'
)
parser.add_argument(
    '--in', nargs='+',
    help='input assay files'
)
parser.add_argument(
    '--reference', nargs='+',
    help='reference assay files'
)
parser.add_argument(
    '--out',
    help='output assay file'
)
parser.add_argument(
    '--out_anchor',
    help='output anchor assay file'
)

multirun_group = parser.add_mutually_exclusive_group(required=False)
multirun_group.add_argument(
    '--multiple_runs', 
    dest='multiple_runs', action='store_true',
    help='assume input files are from different runs (default: %(default)s)'
)
multirun_group.add_argument(
    '--same_run', 
    dest='multiple_runs', action='store_false',
    help='assume input files are from same run (default: True)'
)
parser.set_defaults(multiple_runs=False)

parser.add_argument(
    '--model', choices=['interpolate', 'linear'], default='interpolate',
    help='retention time mapping model (default: %(default)s)'
)

smooth_parameter_group = parser.add_argument_group('smoothing parameters') 
smooth_parameter_group.add_argument(
    '--smooth', choices=['lowess', 'savgol', 'none'], default='lowess',
    help='smoothing method (default: %(default)s)'
)
smooth_parameter_group.add_argument(
    '--lowess_frac', default=0.667, type=float,
    help='LOWESS fraction (default: %(default)s)'
)
smooth_parameter_group.add_argument(
    '--lowess_it', default=0, type=int,
    help='the number of residual-based reweightings in LOWESS (default: %(default)s)'
)
smooth_parameter_group.add_argument(
    '--savgol_window_length', default=7, type=int,
    help='Savitzky-Golay window length (default: %(default)s)'
)
smooth_parameter_group.add_argument(
    '--savgol_polyorder', default=1, type=int,
    help='Savitzky-Golay polyorder (default: %(default)s)'
)

args = parser.parse_args()
assay_files = getattr(args, 'in')
reference_assay_files = args.reference
out_file = args.out
out_anchor_file = args.out_anchor
multiple_runs = args.multiple_runs
model = args.model
smooth = args.smooth
smooth_args = {
    k[len(smooth + '_'):]: v
    for k, v in vars(args).items()
    if k.startswith(smooth + '_') and v is not None
}

    
# %%
import logging

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s %(filename)s: [%(levelname)s] %(message)s'
)

# %%   
if globals().get('assay_files', None) is None or \
    len(assay_files) == 0:
    raise ValueError('no assay files')

if globals().get('reference_assay_files', None) is None or \
    len(reference_assay_files) == 0:
    raise ValueError('no reference assay files')
    
# %%
import os

if globals().get('out_file', None) is None:
    out_file = os.path.splitext(assay_files[0])[0]
    if out_file.endswith('.assay'):
        out_file = out_file[:-len('.assay')]
    if len(assay_files) > 1:
        out_file += '_' + str(len(assay_files))
    out_file += '_rtcalibrated.assay.pickle'
    
if globals().get('out_anchor_file', None) is None:
    out_anchor_file = os.path.splitext(reference_assay_files[0])[0]
    if out_anchor_file.endswith('.assay'):
        out_anchor_file = out_anchor_file[:-len('.assay')]    
    if len(reference_assay_files) > 1:
        out_anchor_file += '_' + str(len(reference_assay_files))
    out_anchor_file += '_rtanchor.assay.pickle'

# %%
logging.info('use model: ' + str(model))

logging.info(
    'use smoothing: ' + str(smooth) + '\n' + \
    '\n'.join((
        k + '=' + str(v) 
        for k, v in smooth_args.items()
        if v is not None
    ))
)
    
# %%
from util import load_pickle, save_pickle
from assay.rtcalibration import RetentionTimeCalibrator

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
reference_assays = []
for reference_assay_file in reference_assay_files:
    logging.info('loading assays: ' + assay_file)  
    
    assay_data = load_pickle(reference_assay_file)
    reference_assays.extend(assay_data)
    
    logging.info('assays loaded: {0}, {1} spectra' \
        .format(reference_assay_file, len(assay_data)))

logging.info('reference assays loaded: {0} spectra totally' \
    .format(len(reference_assays))) 

# %% 
logging.info('calibrating retention time')

calibrator = RetentionTimeCalibrator(
    model=model, 
    smooth=smooth,
    smooth_args=smooth_args
)
calibrator.load_reference(reference_assays)
rt_data = calibrator.calibrate_rt(
    assays, 
    multiple_runs=multiple_runs, 
    inplace=True, 
    return_data=True
)

logging.info('retention time calibrated')

# %%
logging.info('saving assays: {0}' \
    .format(out_file))

save_pickle(assays, out_file)
    
logging.info('assays saved: {0}, {1} spectra' \
    .format(out_file, len(assays)))

# %%
logging.info('saving anchor assays: {0}' \
    .format(out_anchor_file))

anchor_assays = [
    reference_assays[i]
    for i in rt_data['index_reference'] \
        .drop_duplicates().sort_values()
    if i >= 0
]
save_pickle(anchor_assays, out_anchor_file)
    
logging.info('anchor assays saved: {0}, {1} spectra' \
    .format(out_anchor_file, len(anchor_assays)))

# %%
out_rt_file = os.path.splitext(out_file)[0]
if out_rt_file.endswith('.assay'):
    out_rt_file = out_rt_file[:-len('.assay')]
out_rt_file += '.csv' 

logging.info('saving retention time table: {0}' \
             .format(out_rt_file))

rt_data.to_csv(out_rt_file, index=False)

logging.info('retention time saved: {0}' \
             .format(out_rt_file))    

# %%
try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def plot_rt_calibration(rt_data, pdf_path):
    if plt is None:
        raise ImportError("Error: The matplotlib package is required to create a report.")
        
    with PdfPages(pdf_path) as pdf:
        groups = rt_data.groupby(['run'])
              
        if len(groups) > 1:            
            for idx, (run, data) in enumerate(groups):
                if idx % 6 == 0:
                    plt.figure(figsize=(10, 15))
                    plt.subplots_adjust(hspace=.75)
                
                plt.subplot(321 + idx % 6)                            
                        
                plt.scatter(data.rt_old, data.rt_new, marker='.')
                plt.scatter(data.rt_old, data.rt_reference, color='red', marker='D')
                plt.xlabel('Raw RT')
                plt.ylabel('Calibrated RT')
                plt.title(run)
                
                if idx % 6 == 5 or idx == len(groups) - 1:            
                    pdf.savefig()
                    plt.close()
                    
        else:
            plt.figure(figsize=(10, 10))
            
            plt.scatter(rt_data.rt_old, rt_data.rt_new, marker='.')
            plt.scatter(rt_data.rt_old, rt_data.rt_reference, color='red', marker='D')
            plt.xlabel('Raw RT')
            plt.ylabel('Calibrated RT')
            plt.title(rt_data['run'][0])
                
            pdf.savefig()
            plt.close()
         
        
plot_rt_calibration(rt_data, out_rt_file[:-len('.csv')] + '.pdf')


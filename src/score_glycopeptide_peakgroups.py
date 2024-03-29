from pyprophet.data_handling import transform_pi0_lambda, \
    transform_threads
from glycoprophet.runner import GlycoPeptideProphetLearner, \
    GlycoPeptideProphetWeightApplier

import click
from hyperopt import hp


@click.command()
# File handling
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Input file.')
@click.option('--out', 'outfile', type=click.Path(exists=False), help='Output file.')
# Semi-supervised learning
@click.option('--classifier', default='LDA', show_default=True, type=click.Choice(['LDA', 'XGBoost']), help='Either a "LDA" or "XGBoost" classifier is used for semi-supervised learning.')
@click.option('--xgb_autotune/--no-xgb_autotune', default=False, show_default=True, help='XGBoost: Autotune hyperparameters.')
@click.option('--apply_weights', type=click.Path(exists=True), help='Apply score weights file instead of semi-supervised learning.')
@click.option('--xeval_fraction', default=0.5, show_default=True, type=float, help='Data fraction used for cross-validation of semi-supervised learning step.')
@click.option('--xeval_num_iter', default=10, show_default=True, type=int, help='Number of iterations for cross-validation of semi-supervised learning step.')
@click.option('--ss_initial_fdr', default=0.15, show_default=True, type=float, help='Initial FDR cutoff for best scoring targets.')
@click.option('--ss_iteration_fdr', default=0.05, show_default=True, type=float, help='Iteration FDR cutoff for best scoring targets.')
@click.option('--ss_num_iter', default=10, show_default=True, type=int, help='Number of iterations for semi-supervised learning step.')
@click.option('--ss_main_score', default="var_xcorr_shape", show_default=True, type=str, help='Main score to start semi-supervised-learning.')
# Statistics
@click.option('--group_id', default="group_id", show_default=True, type=str, help='Group identifier for calculation of statistics.')
@click.option('--density_estimator', default='gmm', show_default=True, type=click.Choice(['kde', 'gmm']), help='Either kernel density estimation ("kde") or Gaussian mixture model ("gmm") is used for score density estimation.')
@click.option('--grid_size', default=256, show_default=True, type=int, help='Number of d-score cutoffs to build grid coordinates for local FDR calculation.')
@click.option('--parametric/--no-parametric', default=False, show_default=True, help='Do parametric estimation of p-values.')
@click.option('--pfdr/--no-pfdr', default=False, show_default=True, help='Compute positive false discovery rate (pFDR) instead of FDR.')
@click.option('--pi0_lambda', default=[0.1,0.5,0.05], show_default=True, type=(float, float, float), help='Use non-parametric estimation of p-values. Either use <START END STEPS>, e.g. 0.1, 1.0, 0.1 or set to fixed value, e.g. 0.4, 0, 0.', callback=transform_pi0_lambda)
@click.option('--pi0_method', default='bootstrap', show_default=True, type=click.Choice(['smoother', 'bootstrap']), help='Either "smoother" or "bootstrap"; the method for automatically choosing tuning parameter in the estimation of pi_0, the proportion of true null hypotheses.')
@click.option('--pi0_smooth_df', default=3, show_default=True, type=int, help='Number of degrees-of-freedom to use when estimating pi_0 with a smoother.')
@click.option('--pi0_smooth_log_pi0/--no-pi0_smooth_log_pi0', default=False, show_default=True, help='If True and pi0_method = "smoother", pi0 will be estimated by applying a smoother to a scatterplot of log(pi0) estimates against the tuning parameter lambda.')
@click.option('--lfdr_truncate/--no-lfdr_truncate', show_default=True, default=True, help='If True, local FDR values >1 are set to 1.')
@click.option('--lfdr_monotone/--no-lfdr_monotone', show_default=True, default=True, help='If True, local FDR values are non-increasing with increasing d-scores.')
# OpenSWATH options
@click.option('--level', default='ms1ms2', show_default=True, type=click.Choice(['ms2', 'ms1ms2']), help='Either "ms2" or "ms1ms2"; the data level selected for scoring. "ms1ms2 integrates both MS1- and MS2-level scores and can be used instead of "ms2"-level results."')
# TRIC
@click.option('--tric_chromprob/--no-tric_chromprob', default=False, show_default=True, help='Whether chromatogram probabilities for TRIC should be computed.')
# Processing
@click.option('--threads', default=1, show_default=True, type=int, help='Number of threads used for semi-supervised learning. -1 means all available CPUs.', callback=transform_threads)
@click.option('--test/--no-test', default=False, show_default=True, help='Run in test mode with fixed seed.')
def score(infile, outfile, classifier,
          xgb_autotune, apply_weights, xeval_fraction, xeval_num_iter,
          ss_initial_fdr, ss_iteration_fdr, ss_num_iter, ss_main_score,
          group_id, density_estimator, grid_size, parametric, pfdr,
          pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0,
          lfdr_truncate, lfdr_monotone, level, tric_chromprob, threads, test):
    """
    Conduct semi-supervised learning and error-rate estimation for MS1 and MS2-level data.
    """

    if outfile is None:
        outfile = infile
    else:
        outfile = outfile

    # Prepare XGBoost-specific parameters
    xgb_hyperparams = {'autotune': xgb_autotune, 'autotune_num_rounds': 10, 'num_boost_round': 100, 'early_stopping_rounds': 10, 'test_size': 0.33}

    xgb_params = {'eta': 0.3, 'gamma': 0, 'max_depth': 6, 'min_child_weight': 1, 'subsample': 1, 'colsample_bytree': 1, 'colsample_bylevel': 1, 'colsample_bynode': 1, 'lambda': 1, 'alpha': 0, 'scale_pos_weight': 1, 'silent': 1, 'objective': 'binary:logitraw', 'nthread': 1, 'eval_metric': 'auc'}

    xgb_params_space = {'eta': hp.uniform('eta', 0.0, 0.3), 'gamma': hp.uniform('gamma', 0.0, 0.5), 'max_depth': hp.quniform('max_depth', 2, 8, 1), 'min_child_weight': hp.quniform('min_child_weight', 1, 5, 1), 'subsample': 1, 'colsample_bytree': 1, 'colsample_bylevel': 1, 'colsample_bynode': 1, 'lambda': hp.uniform('lambda', 0.0, 1.0), 'alpha': hp.uniform('alpha', 0.0, 1.0), 'scale_pos_weight': 1.0, 'silent': 1, 'objective': 'binary:logitraw', 'nthread': 1, 'eval_metric': 'auc'}

    if not apply_weights:
        GlycoPeptideProphetLearner(
            infile, outfile,
            classifier=classifier,
            xgb_hyperparams=xgb_hyperparams, xgb_params=xgb_params,
            xgb_params_space=xgb_params_space,
            xeval_fraction=xeval_fraction,
            xeval_num_iter=xeval_num_iter,
            ss_initial_fdr=ss_initial_fdr,
            ss_iteration_fdr=ss_iteration_fdr,
            ss_num_iter=ss_num_iter,
            ss_main_score=ss_main_score,
            group_id=group_id,
            density_estimator=density_estimator,
            grid_size=grid_size,
            parametric=parametric, pfdr=pfdr,
            pi0_lambda=pi0_lambda, pi0_method=pi0_method,
            pi0_smooth_df=pi0_smooth_df,
            pi0_smooth_log_pi0=pi0_smooth_log_pi0,
            lfdr_truncate=lfdr_truncate,
            lfdr_monotone=lfdr_monotone,
            level=level,
            tric_chromprob=tric_chromprob,
            threads=threads, test=test
        ).run()
    else:
        GlycoPeptideProphetWeightApplier(
            infile, outfile, apply_weights,
            classifier=classifier,
            xgb_hyperparams=xgb_hyperparams, xgb_params=xgb_params,
            xgb_params_space=xgb_params_space,
            xeval_fraction=xeval_fraction,
            xeval_num_iter=xeval_num_iter,
            ss_initial_fdr=ss_initial_fdr,
            ss_iteration_fdr=ss_iteration_fdr,
            ss_num_iter=ss_num_iter,
            ss_main_score=ss_main_score,
            group_id=group_id,
            density_estimator=density_estimator,
            grid_size=grid_size,
            parametric=parametric, pfdr=pfdr,
            pi0_lambda=pi0_lambda, pi0_method=pi0_method,
            pi0_smooth_df=pi0_smooth_df,
            pi0_smooth_log_pi0=pi0_smooth_log_pi0,
            lfdr_truncate=lfdr_truncate,
            lfdr_monotone=lfdr_monotone,
            level=level,
            tric_chromprob=tric_chromprob,
            threads=threads, test=test
        ).run()


if __name__ == '__main__':
    score()
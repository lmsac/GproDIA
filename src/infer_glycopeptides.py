from pyprophet.data_handling import transform_pi0_lambda
from glycoprophet.level_contexts import infer_glycopeptides

import click


@click.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Input file.')
@click.option('--out', 'outfile', type=click.Path(exists=False), help='Output file.')
@click.option('--context', default='run-specific', show_default=True, type=click.Choice(['run-specific', 'experiment-wide', 'global']), help='Context to estimate glycopeptide-level FDR control.')
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
def glycopeptide(infile, outfile, context, density_estimator, grid_size, parametric, pfdr, pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, lfdr_truncate, lfdr_monotone):
    """
    Infer glycopeptides and conduct error-rate estimation in different contexts.
    """
    if outfile is None:
        outfile = infile
    
    infer_glycopeptides(
        infile, outfile, 
        context=context, 
        density_estimator=density_estimator,
        grid_size=grid_size,
        parametric=parametric, pfdr=pfdr, 
        pi0_lambda=pi0_lambda, pi0_method=pi0_method, 
        pi0_smooth_df=pi0_smooth_df, 
        pi0_smooth_log_pi0=pi0_smooth_log_pi0, 
        lfdr_truncate=lfdr_truncate, 
        lfdr_monotone=lfdr_monotone
    )


if __name__ == '__main__':
    glycopeptide()
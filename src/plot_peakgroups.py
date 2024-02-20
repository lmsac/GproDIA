from glycoprophet.plot import export_peakgroup_plots
import click


@click.command()
@click.option('--in', 'osw_file', required=True, type=click.Path(exists=True), help='Input OSW file.')
@click.option('--mzml', 'mzml_file', required=True, type=click.Path(exists=True), help='Input mzML file.')
@click.option('--swath_windows', 'swath_window_file', type=click.Path(exists=True), help='SWATH window file (TSV).')
@click.option('--transition_groups', 'transition_group_file', type=click.Path(exists=True), help='Transition group-level result file (TSV).')
@click.option('--glycoform/--no-glycoform', 'glycoform', default=False, show_default=True, help='Export glycoform results.')
@click.option('--max_peakgroup_rank', default=1, show_default=True, type=int, help='Filter results to candidate peak groups until maximum peak group rank.')
@click.option('--max_rs_peakgroup_qvalue', default=0.05, show_default=True, type=float, help='Filter results to maximum run-specific peak group-level q-value.')
@click.option('--max_glycoform_qvalue', default=0.05, show_default=True, type=float, help='Filter results to maximum glycoform q-value.')
@click.option('--max_transition_pep', default=0.7, show_default=True, type=float, help='Filter results to maximum transition PEP.')
@click.option('--include_decoy/--exclude_decoy', 'include_decoy', default=False, show_default=True, help='Include decoy results.')

@click.option('--tolerance', default=20, show_default=True, type=float, help='m/z extraction window.')
@click.option('--tolerance_unit', default='ppm', show_default=True, type=click.Choice(['ppm', 'Da']), help='m/z extraction window unit.')
def export_plots(osw_file, mzml_file, swath_window_file, transition_group_file,
                 glycoform, max_peakgroup_rank,
                 max_rs_peakgroup_qvalue, max_glycoform_qvalue, max_transition_pep, include_decoy,
                 tolerance, tolerance_unit):
    """
    Export peak group plots
    """
    export_peakgroup_plots(
        osw_file=osw_file,
        mzml_file=mzml_file,
        swath_window_file=swath_window_file,
        transition_group_file=transition_group_file,
        glycoform=glycoform,
        max_peakgroup_rank=max_peakgroup_rank,
        max_rs_peakgroup_qvalue=max_rs_peakgroup_qvalue,
        max_glycoform_qvalue=max_glycoform_qvalue,
        max_transition_pep=max_transition_pep,
        include_decoy=include_decoy,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit
    )


if __name__ == '__main__':
    export_plots()

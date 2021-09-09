import click

from glycoprophet.glycoform import infer_glycoforms


@click.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Input file.')
@click.option('--out', 'outfile', type=click.Path(exists=False), help='Output file.')
@click.option('--ms1_precursor_scoring/--no-ms1_precursor_scoring', default=True, show_default=True, help='Use MS1 precursor data for glycoform inference.')
@click.option('--ms2_precursor_scoring/--no-ms2_precursor_scoring', default=True, show_default=True, help='Use MS2 precursor data for glycoform inference.')
@click.option('--grouped_fdr/--no-grouped_fdr', default=False, show_default=True, help='[Experimental] Compute grouped FDR instead of pooled FDR to better support data where peak groups are evaluated to originate from very heterogeneous numbers of glycoforms.')
@click.option('--max_precursor_pep', default=1, show_default=True, type=float, help='Maximum PEP to consider scored precursors.')
@click.option('--max_peakgroup_pep', default=0.7, show_default=True, type=float, help='Maximum PEP to consider scored peak groups.')
@click.option('--max_precursor_peakgroup_pep', default=1, show_default=True, type=float, help='Maximum BHM layer 1 integrated precursor peakgroup PEP to consider.')
@click.option('--max_transition_pep', default=0.6, show_default=True, type=float, help='Maximum PEP to consider scored transitions.')
@click.option('--use_glycan_composition/--use_glycan_struct', 'use_glycan_composition', default=True, show_default=True, help='Compute glycoform-level FDR based on glycan composition or struct.')
@click.option('--ms1_mz_window', default=10, show_default=True, type=float, help='MS1 m/z window in Thomson or ppm.')
@click.option('--ms1_mz_window_unit', default='ppm', show_default=True, type=click.Choice(['ppm', 'Da', 'Th']), help='MS1 m/z window unit.')

def glycoform(infile, outfile, 
              ms1_precursor_scoring, ms2_precursor_scoring,
              grouped_fdr,
              max_precursor_pep, max_peakgroup_pep,
              max_precursor_peakgroup_pep,
              max_transition_pep,
              use_glycan_composition,
              ms1_mz_window,
              ms1_mz_window_unit):
    """
    Infer glycoforms after scoring of MS1, MS2 and transition-level data.
    """
    
    if outfile is None:
        outfile = infile
        
    infer_glycoforms(
        infile=infile, outfile=outfile, 
        ms1_precursor_scoring=ms1_precursor_scoring,
        ms2_precursor_scoring=ms2_precursor_scoring,
        grouped_fdr=grouped_fdr,
        max_precursor_pep=max_precursor_pep,
        max_peakgroup_pep=max_peakgroup_pep,
        max_precursor_peakgroup_pep=max_precursor_peakgroup_pep,
        max_transition_pep=max_transition_pep,
        use_glycan_composition=use_glycan_composition,
        ms1_mz_window=ms1_mz_window,
        ms1_mz_window_unit=ms1_mz_window_unit
    )


if __name__ == '__main__':
    glycoform()

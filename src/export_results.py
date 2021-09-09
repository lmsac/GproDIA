from glycoprophet.export import export_tsv, export_score_plots
import click


@click.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Input file.')
@click.option('--out', 'outfile', type=click.Path(exists=False), help='Output TSV/CSV (matrix, legacy_split, legacy_merged) file.')
@click.option('--format', default='legacy_split', show_default=True, type=click.Choice(['matrix', 'legacy_split', 'legacy_merged','score_plots']), help='Export format, either matrix, legacy_split/legacy_merged (mProphet/PyProphet) or score_plots format.')
@click.option('--csv/--no-csv', 'outcsv', default=False, show_default=True, help='Export CSV instead of TSV file.')
@click.option('--transition_quantification/--no-transition_quantification', default=True, show_default=True, help='[format: legacy] Report aggregated transition-level quantification.')
@click.option('--max_transition_pep', default=0.7, show_default=True, type=float, help='[format: legacy] Maximum PEP to retain scored transitions for quantification (requires transition-level scoring).')

@click.option('--glycoform/--no-glycoform', 'glycoform', default=False, show_default=True, help='[format: matrix/legacy] Export glycoform results.')
@click.option('--glycoform_match_precursor', default='glycan_composition', show_default=True, type=click.Choice(['exact', 'glycan_composition', 'none']), help='[format: matrix/legacy] Export glycoform results with glycan matched with precursor-level results.')
@click.option('--max_glycoform_pep', default=1, show_default=True, type=float, help='[format: matrix/legacy] Filter results to maximum glycoform PEP.')
@click.option('--max_glycoform_qvalue', default=0.05, show_default=True, type=float, help='[format: matrix/legacy] Filter results to maximum glycoform q-value.')

@click.option('--max_rs_peakgroup_qvalue', default=0.05, show_default=True, type=float, help='[format: matrix/legacy] Filter results to maximum run-specific peak group-level q-value.')
@click.option('--glycopeptide/--no-glycopeptide', default=True, show_default=True, help='Append glycopeptide-level error-rate estimates if available.')
@click.option('--max_global_glycopeptide_qvalue', default=0.01, show_default=True, type=float, help='[format: matrix/legacy] Filter results to maximum global glycopeptide-level q-value.')
def export(infile, outfile, format, outcsv, 
           transition_quantification, max_transition_pep, 
           glycoform, glycoform_match_precursor,
           max_glycoform_pep, max_glycoform_qvalue,
           max_rs_peakgroup_qvalue, 
           glycopeptide, max_global_glycopeptide_qvalue):
    """
    Export TSV/CSV tables
    """
    if format == "score_plots":
        export_score_plots(infile)
    else:
        if outfile is None:
            if outcsv:
                outfile = infile.split(".osw")[0] + ".csv"
            else:
                outfile = infile.split(".osw")[0] + ".tsv"
        else:
            outfile = outfile

        export_tsv(
            infile, outfile, 
            format=format, outcsv=outcsv, 
            transition_quantification=transition_quantification, 
            max_transition_pep=max_transition_pep, 
            glycoform=glycoform, 
            glycoform_match_precursor=glycoform_match_precursor,
            max_glycoform_pep=max_glycoform_pep, 
            max_glycoform_qvalue=max_glycoform_qvalue,
            max_rs_peakgroup_qvalue=max_rs_peakgroup_qvalue, 
            glycopeptide=glycopeptide,
            max_global_glycopeptide_qvalue=max_global_glycopeptide_qvalue,
        )

    
if __name__ == '__main__':
    export()
from glycoprophet.level_contexts import backpropagate_scores
import click


@click.command()
@click.option('--in', 'infile', required=True, type=click.Path(exists=True), help='Single run input file.')
@click.option('--out','outfile', type=click.Path(exists=False), help='Single run (with multi-run scores) output file.')
@click.option('--apply_scores', required=True, type=click.Path(exists=True), help='Multi-run scores file to apply.')
def backpropagate(infile, outfile, apply_scores):
    """
    Backpropagate multi-run glycopeptide scores to single files
    """

    if outfile is None:
        outfile = infile
    else:
        outfile = outfile

    backpropagate_scores(infile, outfile, apply_scores)


if __name__ == '__main__':
    backpropagate()
from glycoprophet.oxonium import validate_oxonium_ions
import click


@click.command()
@click.option('--in', 'report_file', required=True, type=click.Path(exists=True), help='Input TSV/CSV file.')
@click.option('--mzml', 'mzml_files', nargs=0, type=click.Path(exists=True), help='Input mzML files.')
@click.argument('mzml_files', required=True, nargs=-1, type=click.Path(exists=True))
@click.option('--out', 'out_file', type=click.Path(exists=False), help='Output TSV/CSV file.')
@click.option('--mz_tolerance', default=20, show_default=True, type=float, help='MS2 m/z window in Thomson or ppm.')
@click.option('--mz_tolerance_unit', default='ppm', show_default=True, type=click.Choice(['ppm', 'Da', 'Th']), help='MS2 m/z window unit.')
@click.option('--remove_background/--no-remove_background', 'background_estimator', default=True, show_default=True, help='Remove background signals.')
@click.option('--absolute_intensity', default=10, show_default=True, type=float, help='Minimum absolute intensity to consider a peak.')
@click.option('--relative_intensity', default=0.05, show_default=True, type=float, help='Minimum relative intensity to consider a peak.')
def oxonium(report_file, mzml_files, out_file,
            mz_tolerance, mz_tolerance_unit,
            background_estimator, absolute_intensity, relative_intensity):    
    validate_oxonium_ions(
        report_file, mzml_files, out_file, 
        mz_tolerance=mz_tolerance, mz_tolerance_unit=mz_tolerance_unit,
        background_estimator=background_estimator,
        absolute_intensity=absolute_intensity,
        relative_intensity=relative_intensity
    )
    

if __name__ == '__main__':
    oxonium()
    
    
from pyprophet.data_handling import check_sqlite_table
import pandas as pd
import sqlite3
import click
from shutil import copyfile

from .stats import ErrorStatisticsCalculator
from .report import save_report

def statistics_report(data, outfile, context, analyte, 
                      **kwargs):

    error_stat = ErrorStatisticsCalculator(data, **kwargs)    
    result, pi0 = error_stat.error_statistics()
    
    # print summary table
    summary = result.get('summary_statistics', None)
    if summary is not None:
        summary = summary.get('total', None)
    
    if summary is not None:
        click.echo("=" * 80)
        click.echo(summary)
        click.echo("=" * 80)

    if context == 'run-specific':
        outfile = outfile + "_" + str(data['run_id'].unique()[0])

    # export PDF report
    save_report(
        outfile + "_" + context + "_" + analyte + ".pdf", 
        outfile + ": \n" + context + " " + analyte + "-level error-rate control", 
        result['scored_table'],
        result['final_statistics'],
        pi0
    )
    return result['scored_table']
    

def infer_glycopeptides(infile, outfile, context,
                        density_estimator,
                        grid_size,
                        parametric, pfdr,
                        pi0_lambda,
                        pi0_method, pi0_smooth_df, 
                        pi0_smooth_log_pi0,
                        lfdr_truncate, lfdr_monotone, 
                        # lfdr_transformation, lfdr_adj, lfdr_eps
                        ):
    con = sqlite3.connect(infile)
    
    if not check_sqlite_table(con, "SCORE_MS2") or \
        not check_sqlite_table(con, "SCORE_MS2_PART_PEPTIDE") or \
        not check_sqlite_table(con, "SCORE_MS2_PART_GLYCAN"):
        raise click.ClickException("Apply scoring to MS2-level data before running glycopeptide-level scoring.")
        
    if context not in ['global','experiment-wide','run-specific']:
        raise click.ClickException("Unspecified context selected.")
        
    if context == 'global':
        run_id = 'NULL'
        group_id = 'GLYCOPEPTIDE.ID'
    else:
        run_id = 'RUN_ID'
        group_id = 'RUN_ID || "_" || GLYCOPEPTIDE.ID'
        
    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_glycopeptide_glycopeptide_id ON GLYCOPEPTIDE (ID);
CREATE INDEX IF NOT EXISTS idx_precursor_glycopeptide_mapping_glycopeptide_id ON PRECURSOR_GLYCOPEPTIDE_MAPPING (GLYCOPEPTIDE_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_glycopeptide_mapping_precursor_id ON PRECURSOR_GLYCOPEPTIDE_MAPPING (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_part_peptide_feature_id ON SCORE_MS2_PART_PEPTIDE (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_part_glycan_feature_id ON SCORE_MS2_PART_GLYCAN (FEATURE_ID);
''')
    
    data = pd.read_sql_query('''
SELECT %s AS RUN_ID,
       %s AS GROUP_ID,
       GLYCOPEPTIDE.ID AS GLYCOPEPTIDE_ID,
       GLYCOPEPTIDE.DECOY_PEPTIDE,
       GLYCOPEPTIDE.DECOY_GLYCAN,
       SCORE_MS2.SCORE AS d_score_combined,
       SCORE_MS2_PART_PEPTIDE.SCORE AS d_score_peptide,
       SCORE_MS2_PART_GLYCAN.SCORE AS d_score_glycan,
       "%s" AS CONTEXT
FROM GLYCOPEPTIDE
INNER JOIN PRECURSOR_GLYCOPEPTIDE_MAPPING ON GLYCOPEPTIDE.ID = PRECURSOR_GLYCOPEPTIDE_MAPPING.GLYCOPEPTIDE_ID
INNER JOIN PRECURSOR ON PRECURSOR_GLYCOPEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
INNER JOIN SCORE_MS2_PART_PEPTIDE ON FEATURE.ID = SCORE_MS2_PART_PEPTIDE.FEATURE_ID
INNER JOIN SCORE_MS2_PART_GLYCAN ON FEATURE.ID = SCORE_MS2_PART_GLYCAN.FEATURE_ID
GROUP BY GROUP_ID
HAVING MAX(d_score_combined)
ORDER BY d_score_combined DESC
''' % (run_id, group_id, context), con)
    
    data.columns = [col.lower() for col in data.columns]
    
    if context == 'run-specific':
        data = data.groupby('run_id').apply(
            statistics_report, 
            outfile, context, 'glycopeptide',
            density_estimator=density_estimator,
            grid_size=grid_size,
            parametric=parametric, pfdr=pfdr, 
            pi0_lambda=pi0_lambda, pi0_method=pi0_method, 
            pi0_smooth_df=pi0_smooth_df, 
            pi0_smooth_log_pi0=pi0_smooth_log_pi0, 
            lfdr_truncate=lfdr_truncate, 
            lfdr_monotone=lfdr_monotone, 
            # lfdr_transformation=lfdr_transformation, 
            # lfdr_adj=lfdr_adj, lfdr_eps=lfdr_eps
        ).reset_index()

    elif context in ['global', 'experiment-wide']:
        data = statistics_report(
            data, outfile, context, 'glycopeptide',
            density_estimator=density_estimator,
            grid_size=grid_size,
            parametric=parametric, pfdr=pfdr, 
            pi0_lambda=pi0_lambda, pi0_method=pi0_method, 
            pi0_smooth_df=pi0_smooth_df, 
            pi0_smooth_log_pi0=pi0_smooth_log_pi0, 
            lfdr_truncate=lfdr_truncate, 
            lfdr_monotone=lfdr_monotone, 
            # lfdr_transformation=lfdr_transformation, 
            # lfdr_adj=lfdr_adj, lfdr_eps=lfdr_eps
        )
    
    if infile != outfile:
        copyfile(infile, outfile)
    
    con = sqlite3.connect(outfile)

    c = con.cursor()
    c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="SCORE_GLYCOPEPTIDE"')
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_GLYCOPEPTIDE WHERE CONTEXT =="%s"' % context)
    c.fetchall()
    c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="SCORE_GLYCOPEPTIDE_PART_PEPTIDE"')
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_GLYCOPEPTIDE_PART_PEPTIDE WHERE CONTEXT =="%s"' % context)
    c.fetchall()
    c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="SCORE_GLYCOPEPTIDE_PART_GLYCAN"')
    if c.fetchone()[0] == 1:
        c.execute('DELETE FROM SCORE_GLYCOPEPTIDE_PART_GLYCAN WHERE CONTEXT =="%s"' % context)
    c.fetchall()

    df = data[['context','run_id','glycopeptide_id','d_score_combined','q_value','pep']]
    df.columns = ['CONTEXT','RUN_ID','GLYCOPEPTIDE_ID','SCORE','QVALUE','PEP']
    table = "SCORE_GLYCOPEPTIDE"
    df.to_sql(table, con, index=False, dtype={"RUN_ID": "INTEGER"}, if_exists='append')
    
    for part in ['peptide', 'glycan']:
        df = data[['context','run_id','glycopeptide_id','d_score_' + part,'pep_' + part]]
        df.columns = ['CONTEXT','RUN_ID','GLYCOPEPTIDE_ID','SCORE','PEP']
        table = "SCORE_GLYCOPEPTIDE_PART_" + part.upper()
        df.to_sql(table, con, index=False, dtype={"RUN_ID": "INTEGER"}, if_exists='append')

    con.close()
    

def backpropagate_scores(infile, outfile, apply_scores):
    if infile != outfile:
        copyfile(infile, outfile)

    score_con = sqlite3.connect(apply_scores)
    glycopeptide_present = \
        check_sqlite_table(score_con, "SCORE_GLYCOPEPTIDE") and \
        check_sqlite_table(score_con, "SCORE_GLYCOPEPTIDE_PART_PEPTIDE") and \
        check_sqlite_table(score_con, "SCORE_GLYCOPEPTIDE_PART_GLYCAN")
    score_con.close()
    if not glycopeptide_present:
        raise click.ClickException('Backpropagation requires glycopeptide-level contexts.')

    script = list()
    script.append('PRAGMA synchronous = OFF;')
    script.append('DROP TABLE IF EXISTS SCORE_GLYCOPEPTIDE;')
    script.append('DROP TABLE IF EXISTS SCORE_GLYCOPEPTIDE_PART_PEPTIDE;')
    script.append('DROP TABLE IF EXISTS SCORE_GLYCOPEPTIDE_PART_GLYCAN;')

    if glycopeptide_present:
        script.append('CREATE TABLE SCORE_GLYCOPEPTIDE (CONTEXT TEXT, RUN_ID INTEGER, GLYCOPEPTIDE_ID INTEGER, SCORE REAL, QVALUE REAL, PEP REAL);')
        script.append('CREATE TABLE SCORE_GLYCOPEPTIDE_PART_PEPTIDE (CONTEXT TEXT, RUN_ID INTEGER, GLYCOPEPTIDE_ID INTEGER, SCORE REAL, PEP REAL);')
        script.append('CREATE TABLE SCORE_GLYCOPEPTIDE_PART_GLYCAN (CONTEXT TEXT, RUN_ID INTEGER, GLYCOPEPTIDE_ID INTEGER, SCORE REAL, PEP REAL);')
        
    script.append('ATTACH DATABASE "{}" AS sdb;'.format(apply_scores))
    insert_table_fmt = 'INSERT INTO {0}\nSELECT *\nFROM sdb.{0};'
    if glycopeptide_present:
        script.append(insert_table_fmt.format('SCORE_GLYCOPEPTIDE'))
        script.append(insert_table_fmt.format('SCORE_GLYCOPEPTIDE_PART_PEPTIDE'))
        script.append(insert_table_fmt.format('SCORE_GLYCOPEPTIDE_PART_GLYCAN'))
        
    conn = sqlite3.connect(outfile)
    c = conn.cursor()
    c.executescript('\n'.join(script))
    conn.commit()
    conn.close()

    click.echo("Info: All multi-run data was backpropagated.")
    
    
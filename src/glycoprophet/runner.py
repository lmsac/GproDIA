import click
import os
import time
import warnings
import pandas as pd
import numpy as np
import sqlite3

from pyprophet.pyprophet import PyProphet
from pyprophet.data_handling import check_sqlite_table
from shutil import copyfile

from .runnerbase import ProphetRunnerBase
from .scoring import LDAGlycoPeptideScorer
from .stats import ErrorStatisticsCalculator
from .report import save_report


def read_osw_ms2_feature(infile, include_ms1=False):
    con = sqlite3.connect(infile)

    if not check_sqlite_table(con, "FEATURE_MS2"):
        raise click.ClickException("MS2-level feature table not present in file.")

    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_feature_ms2_feature_id ON FEATURE_MS2 (FEATURE_ID);
''')

    table = pd.read_sql_query('''
SELECT *,
       RUN_ID || '_' || PRECURSOR_ID AS GROUP_ID
FROM FEATURE_MS2
INNER JOIN
  (SELECT RUN_ID,
          ID,
          PRECURSOR_ID,
          EXP_RT
   FROM FEATURE) AS FEATURE ON FEATURE_ID = FEATURE.ID
INNER JOIN
  (SELECT ID,
          CHARGE AS PRECURSOR_CHARGE,
          DECOY
   FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID

INNER JOIN
  (SELECT PRECURSOR_ID AS ID,
          DECOY_PEPTIDE,
          DECOY_GLYCAN
   FROM PRECURSOR_GLYCOPEPTIDE_MAPPING
   INNER JOIN GLYCOPEPTIDE
   ON PRECURSOR_GLYCOPEPTIDE_MAPPING.GLYCOPEPTIDE_ID == GLYCOPEPTIDE.ID) AS DECOY
  ON FEATURE.PRECURSOR_ID = DECOY.ID

INNER JOIN
  (SELECT PRECURSOR_ID AS ID,
          COUNT(*) AS TRANSITION_COUNT
   FROM TRANSITION_PRECURSOR_MAPPING
   INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID
   WHERE DETECTING==1
   GROUP BY PRECURSOR_ID) AS VAR_TRANSITION_SCORE ON FEATURE.PRECURSOR_ID = VAR_TRANSITION_SCORE.ID
ORDER BY RUN_ID,
         PRECURSOR.ID ASC,
         FEATURE.EXP_RT ASC;
''', con)

    if include_ms1:
        if not check_sqlite_table(con, "FEATURE_MS1"):
            raise click.ClickException("MS1-level feature table not present in file.")
        ms1_table = pd.read_sql_query('SELECT * FROM FEATURE_MS1;', con)

        ms1_scores = [c for c in ms1_table.columns if c.startswith("VAR_")]
        ms1_table = ms1_table[['FEATURE_ID'] + ms1_scores]
        ms1_table.columns = ['FEATURE_ID'] + ["VAR_MS1_" + s.split("VAR_")[1] for s in ms1_scores]

        table = pd.merge(table, ms1_table, how='left', on='FEATURE_ID')

    table.columns = [col.lower() for col in table.columns]

    return table


def read_osw_ms1_feature(infile):
    con = sqlite3.connect(infile)

    if not check_sqlite_table(con, "FEATURE_MS1"):
        raise click.ClickException("MS1-level feature table not present in file.")

    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_feature_ms1_feature_id ON FEATURE_MS1 (FEATURE_ID);
''')

    table = pd.read_sql_query('''
SELECT *,
   RUN_ID || '_' || PRECURSOR_ID AS GROUP_ID
FROM FEATURE_MS1
INNER JOIN
  (SELECT RUN_ID,
      ID,
      PRECURSOR_ID,
      EXP_RT
   FROM FEATURE) AS FEATURE ON FEATURE_ID = FEATURE.ID
INNER JOIN
  (SELECT ID,
      CHARGE AS PRECURSOR_CHARGE,
      DECOY
   FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID

INNER JOIN
      (SELECT PRECURSOR_ID AS ID,
              DECOY_PEPTIDE,
              DECOY_GLYCAN
       FROM PRECURSOR_GLYCOPEPTIDE_MAPPING
       INNER JOIN GLYCOPEPTIDE
       ON PRECURSOR_GLYCOPEPTIDE_MAPPING.GLYCOPEPTIDE_ID == GLYCOPEPTIDE.ID) AS DECOY
      ON FEATURE.PRECURSOR_ID = DECOY.ID

ORDER BY RUN_ID,
      PRECURSOR.ID ASC,
      FEATURE.EXP_RT ASC;
''', con)

    table.columns = [col.lower() for col in table.columns]

    return table


class GlycoPeptideProphetRunner(ProphetRunnerBase):
    def __init__(self, infile, outfile, classifier,
                 xgb_hyperparams, xgb_params, xgb_params_space,
                 xeval_fraction, xeval_num_iter,
                 ss_initial_fdr, ss_iteration_fdr, ss_num_iter, ss_main_score,
                 group_id, density_estimator, grid_size,
                 parametric, pfdr,
                 pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0,
                 lfdr_truncate, lfdr_monotone, # lfdr_transformation,
                 # lfdr_adj, lfdr_eps,
                 level, tric_chromprob, threads, test):
        super(GlycoPeptideProphetRunner, self).__init__(
            infile=infile,
            outfile=outfile,
            classifier=classifier,
            ss_main_score=ss_main_score,
            level=level
        )

        self.xgb_hyperparams = xgb_hyperparams
        self.xgb_params = xgb_params
        self.xgb_params_space = xgb_params_space
        self.xeval_fraction = xeval_fraction
        self.xeval_num_iter = xeval_num_iter
        self.ss_initial_fdr = ss_initial_fdr
        self.ss_iteration_fdr = ss_iteration_fdr
        self.ss_num_iter = ss_num_iter

        self.group_id = group_id
        self.density_estimator = density_estimator
        self.grid_size = grid_size
        self.parametric = parametric
        self.pfdr = pfdr
        self.pi0_lambda = pi0_lambda
        self.pi0_method = pi0_method
        self.pi0_smooth_df = pi0_smooth_df
        self.pi0_smooth_log_pi0 = pi0_smooth_log_pi0
        self.lfdr_truncate = lfdr_truncate
        self.lfdr_monotone = lfdr_monotone
        self.lfdr_transformation = 'probit' # lfdr_transformation
        self.lfdr_adj = 1.5 # lfdr_adj
        self.lfdr_eps = np.power(10.0, -8) # lfdr_eps

        self.tric_chromprob = tric_chromprob
        self.threads = threads
        self.test = test


    def load_osw(self):
        if self.level == "ms2" or self.level == "ms1ms2":
            table = read_osw_ms2_feature(
                self.infile,
                include_ms1=self.level == "ms1ms2"
            )
        elif self.level == "ms1":
            table = read_osw_ms1_feature(self.infile)
        else:
            raise click.ClickException("Unspecified data level selected.")
        return table


    def score(self):
        start_at_peptide = time.time()
        click.echo("Info: scoring peptide part")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            (result_peptide, scorer_peptide, weights_peptide) = \
                self.partial_score(part='peptide')

        end_at_peptide = time.time()
        needed = end_at_peptide - start_at_peptide
        seconds = int(needed)
        msecs = int(1000 * (needed - seconds))
        click.echo("Info: peptide part scored: %d seconds and %d msecs" % (seconds, msecs))

        start_at_glycan = time.time()
        click.echo("Info: scoring glycan part")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            (result_glycan, scorer_glycan, weights_glycan) = \
                self.partial_score(part='glycan')

        end_at_glycan = time.time()
        needed = end_at_glycan - start_at_glycan
        seconds = int(needed)
        msecs = int(1000 * (needed - seconds))
        click.echo("Info: glycan part scored: %d seconds and %d msecs" % (seconds, msecs))

        start_at_combined = time.time()
        click.echo("Info: calculating combined scores")
        (result_combined, weights_combined) = \
            self.combined_score(result_peptide, result_glycan)

        end_at_combined= time.time()
        if isinstance(weights_combined, pd.DataFrame):
            click.echo(weights_combined)
        needed = end_at_combined - start_at_combined
        seconds = int(needed)
        msecs = int(1000 * (needed - seconds))
        click.echo("Info: combined scores calculated: %d seconds and %d msecs" % (seconds, msecs))

        start_at_stats = time.time()
        click.echo("Info: calculating error statistics")
        error_stat = ErrorStatisticsCalculator(
            result_combined,
            density_estimator=self.density_estimator,
            grid_size=self.grid_size,
            parametric=self.parametric,
            pfdr=self.pfdr,
            pi0_lambda=self.pi0_lambda, pi0_method=self.pi0_method,
            pi0_smooth_df=self.pi0_smooth_df,
            pi0_smooth_log_pi0=self.pi0_smooth_log_pi0,
            lfdr_truncate=self.lfdr_truncate,
            lfdr_monotone=self.lfdr_monotone,
            lfdr_transformation=self.lfdr_transformation,
            lfdr_adj=self.lfdr_adj, lfdr_eps=self.lfdr_eps,
            tric_chromprob=self.tric_chromprob,
        )
        result, pi0 = error_stat.error_statistics()

        end_at_stats = time.time()
        needed = end_at_stats - start_at_stats
        seconds = int(needed)
        msecs = int(1000 * (needed - seconds))
        click.echo("Info: error statistics finished: %d seconds and %d msecs" % (seconds, msecs))

        self.print_summary(result)

        if all((
            isinstance(w, pd.DataFrame)
            for w in [weights_peptide, weights_glycan, weights_combined]
        )):
            weights = pd.concat((
                weights_peptide.assign(part='peptide'),
                weights_glycan.assign(part='glycan'),
                weights_combined.assign(part='combined')
            ), ignore_index=True)
        else:
            weights = {
               'peptide': weights_peptide,
               'glycan': weights_glycan,
               'combined': weights_combined
            }

        return result, weights, pi0


    def partial_score(self, part):
        raise NotImplementedError('partial_score')

    def combined_score(self, result_peptide, result_glycan):
        raise NotImplementedError('total_score')


    def extra_writes(self):
        d = super(GlycoPeptideProphetRunner, self).extra_writes()
        d.update({
            'report_path': os.path.join(self.prefix + '_report.pdf')
        })
        return d

    def save_tsv_results(self, result, extra_writes, pi0):
        super(GlycoPeptideProphetRunner, self) \
            .save_tsv_results(result, extra_writes, pi0)

        output_path = extra_writes.get("output_path", None)
        report_path = extra_writes.get("report_path", None)
        if output_path is not None:
            output_path= report_path
        if report_path is not None:
            save_report(
                extra_writes.get("report_path"),
                output_path + ': ' + self.level + '-level scoring',
                result['scored_table'],
                result['final_statistics'],
                pi0
            )
            click.echo("Info: %s written." % extra_writes.get("report_path"))


    def save_osw_results(self, result, extra_writes, pi0):
        if self.infile != self.outfile:
            copyfile(self.infile, self.outfile)

        con = sqlite3.connect(self.outfile)

        if self.level == "ms2" or self.level == "ms1ms2":
            c = con.cursor()
            c.execute('DROP TABLE IF EXISTS SCORE_MS2;')
            c.execute('DROP TABLE IF EXISTS SCORE_MS2_PART_PEPTIDE;')
            c.execute('DROP TABLE IF EXISTS SCORE_MS2_PART_GLYCAN;')
            con.commit()
            c.fetchall()

            table = "SCORE_MS2"

        elif self.level == "ms1":
            c = con.cursor()
            c.execute('DROP TABLE IF EXISTS SCORE_MS1;')
            c.execute('DROP TABLE IF EXISTS SCORE_MS1_PART_PEPTIDE;')
            c.execute('DROP TABLE IF EXISTS SCORE_MS1_PART_GLYCAN;')
            con.commit()
            c.fetchall()

            table = "SCORE_MS1"

        df = result['scored_table']
        if 'h_score' in df.columns:
            df = df[['feature_id','d_score_combined','h_score','h0_score','peak_group_rank','q_value','pep']]
            df.columns = ['FEATURE_ID','SCORE','HSCORE','H0SCORE','RANK','QVALUE','PEP']
        else:
            df = df[['feature_id','d_score_combined','peak_group_rank','q_value','pep']]
            df.columns = ['FEATURE_ID','SCORE','RANK','QVALUE','PEP']
        df.to_sql(table, con, index=False)

        for part in ['peptide', 'glycan']:
            df = result['scored_table']
            df = df[['feature_id','d_score_' + part,'pep_' + part]]
            df.columns = ['FEATURE_ID','SCORE','PEP']
            df.to_sql(table + "_PART_" + part.upper(), con, index=False)

        con.close()
        click.echo("Info: %s written." % self.outfile)

        if result.get('final_statistics', None) is not None:
            save_report(
                os.path.join(self.prefix + "_" + self.level + "_report.pdf"),
                self.outfile + ': ' + self.level + '-level scoring',
                result['scored_table'],
                result['final_statistics'],
                pi0
            )
            click.echo("Info: %s written." % os.path.join(self.prefix + "_" + self.level + "_report.pdf"))


class GlycoPeptideProphetLearner(GlycoPeptideProphetRunner):
    def partial_score(self, part):
        if part != 'peptide' and part != 'glycan':
            raise click.ClickException("Unspecified scoring part selected.")

        table = self.table
        if 'decoy' in table.columns:
            table = table.drop(columns=['decoy'])
        table = table.rename(columns={'decoy_' + part: 'decoy'})

        (result, scorer, weights) = PyProphet(
            classifier=self.classifier,
            xgb_hyperparams=self.xgb_hyperparams, xgb_params=self.xgb_params,
            xgb_params_space=self.xgb_params_space,
            xeval_fraction=self.xeval_fraction,
            xeval_num_iter=self.xeval_num_iter,
            ss_initial_fdr=self.ss_initial_fdr,
            ss_iteration_fdr=self.ss_iteration_fdr,
            ss_num_iter=self.ss_num_iter,
            group_id=self.group_id,
            parametric=self.parametric, pfdr=self.pfdr,
            pi0_lambda=self.pi0_lambda, pi0_method=self.pi0_method,
            pi0_smooth_df=self.pi0_smooth_df,
            pi0_smooth_log_pi0=self.pi0_smooth_log_pi0,
            lfdr_truncate=self.lfdr_truncate,
            lfdr_monotone=self.lfdr_monotone,
            lfdr_transformation=self.lfdr_transformation,
            lfdr_adj=self.lfdr_adj, lfdr_eps=self.lfdr_eps,
            tric_chromprob=self.tric_chromprob,
            threads=self.threads, test=self.test,
        ).learn_and_apply(table)

        return (result, scorer, weights)


    def combined_score(self, result_peptide, result_glycan):
        table = pd.merge(
            result_peptide.scored_tables[[
                self.group_id, 'feature_id', 'run_id', 'precursor_id',
                'd_score', 'decoy'
            ]],
            result_glycan.scored_tables[[
                self.group_id, 'feature_id', 'run_id', 'precursor_id',
                'd_score', 'decoy'
            ]],
            on=[self.group_id, 'feature_id', 'run_id', 'precursor_id'],
            suffixes=['_peptide', '_glycan']
        )

        scorer = LDAGlycoPeptideScorer(group_id=self.group_id)
        result = scorer.learn_and_apply(table)
        return (result['scored_table'], result['weights'])


    def extra_writes(self):
        d = super(GlycoPeptideProphetLearner, self).extra_writes()
        d.update({
            'trained_weights_path': os.path.join(self.prefix + '_weights.csv'),
            'trained_model_path_ms1': os.path.join(self.prefix + '_ms1_model.bin'),
            'trained_model_path_ms1ms2': os.path.join(self.prefix + '_ms1ms2_model.bin'),
            'trained_model_path_ms2': os.path.join(self.prefix + '_ms2_model.bin')
        })
        return d


class GlycoPeptideProphetWeightApplier(GlycoPeptideProphetRunner):
    def __init__(self, infile, outfile, apply_weights, **kwargs):
        super(GlycoPeptideProphetWeightApplier, self) \
            .__init__(infile, outfile, **kwargs)

        self.persisted_weights = self.load_weights(apply_weights)


    def partial_score(self, part):
        if part != 'peptide' and part != 'glycan':
            raise click.ClickException("Unspecified scoring part selected.")

        if isinstance(self.persisted_weights, pd.DataFrame):
            weights = self.persisted_weights \
                .loc[self.persisted_weights['part'] == part] \
                .drop(columns=['part'])
        else:
            weights = self.persisted_weights[part]

        table = self.table
        if 'decoy' in table.columns:
            table = table.drop(columns=['decoy'])
        table = table.rename(columns={'decoy_' + part: 'decoy'})

        (result, scorer, weights) = PyProphet(
            classifier=self.classifier,
            xgb_hyperparams=self.xgb_hyperparams, xgb_params=self.xgb_params,
            xgb_params_space=self.xgb_params_space,
            xeval_fraction=self.xeval_fraction,
            xeval_num_iter=self.xeval_num_iter,
            ss_initial_fdr=self.ss_initial_fdr,
            ss_iteration_fdr=self.ss_iteration_fdr,
            ss_num_iter=self.ss_num_iter,
            group_id=self.group_id,
            parametric=self.parametric, pfdr=self.pfdr,
            pi0_lambda=self.pi0_lambda, pi0_method=self.pi0_method,
            pi0_smooth_df=self.pi0_smooth_df,
            pi0_smooth_log_pi0=self.pi0_smooth_log_pi0,
            lfdr_truncate=self.lfdr_truncate,
            lfdr_monotone=self.lfdr_monotone,
            lfdr_transformation=self.lfdr_transformation,
            lfdr_adj=self.lfdr_adj, lfdr_eps=self.lfdr_eps,
            tric_chromprob=self.tric_chromprob,
            threads=self.threads, test=self.test,
        ).apply_weights(table, weights)

        return (result, scorer, weights)


    def combined_score(self, result_peptide, result_glycan):
        if isinstance(self.persisted_weights, pd.DataFrame):
            weights = self.persisted_weights \
                .loc[self.persisted_weights['part'] == 'combined'] \
                .drop(columns=['part'])
        else:
            weights = self.persisted_weights['combined']

        table = pd.merge(
            result_peptide.scored_tables[[
                self.group_id, 'feature_id', 'run_id', 'precursor_id',
                'd_score', 'decoy'
            ]],
            result_glycan.scored_tables[[
                self.group_id, 'feature_id', 'run_id', 'precursor_id',
                'd_score', 'decoy'
            ]],
            on=[self.group_id, 'feature_id', 'run_id', 'precursor_id'],
            suffixes=['_peptide', '_glycan']
        )

        scorer = LDAGlycoPeptideScorer()
        result = scorer.apply_weights(table, weights)
        return (result['scored_table'], result['weights'])


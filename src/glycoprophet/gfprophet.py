import click
import os
import time
import warnings
import pandas as pd
import sqlite3

from pyprophet.pyprophet import PyProphet
from pyprophet.data_handling import check_sqlite_table
from shutil import copyfile

from .runnerbase import ProphetRunnerBase
from .report import save_report_pyprophet as save_report



def read_osw_ms1_feature(infile, max_peakgroup_rank):
    con = sqlite3.connect(infile)
    
    if not check_sqlite_table(con, "SCORE_MS2"):
        raise click.ClickException("MS1-level scoring for glycoform inference requires prior MS2 or MS1MS2-level scoring. Please run 'pyprophet score --level=ms2' or 'pyprophet score --level=ms1ms2' on this file first.")
    if not check_sqlite_table(con, "FEATURE_MS1"):
        raise click.ClickException("MS1-level feature table not present in file.")

    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_feature_ms1_feature_id ON FEATURE_MS1 (FEATURE_ID);
''')

    table = pd.read_sql_query('''
SELECT DECOY.*,
       FEATURE_MS1.*, 
       FEATURE.*,
       PRECURSOR.*,
       RUN_ID || '_' || PRECURSOR_ID AS GROUP_ID
FROM FEATURE_MS1
INNER JOIN
  (SELECT RUN_ID,
      ID,
      PRECURSOR_ID,
      EXP_RT
   FROM FEATURE) AS FEATURE ON FEATURE_MS1.FEATURE_ID = FEATURE.ID
  
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID

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

WHERE RANK <= %s
ORDER BY RUN_ID,
      PRECURSOR.ID ASC,
      FEATURE.EXP_RT ASC;
''' % (max_peakgroup_rank), con)
      
    table.columns = [col.lower() for col in table.columns]
    
    return table 


def read_osw_transition_feature(infile, max_peakgroup_rank, max_peakgroup_pep, max_transition_isotope_overlap, min_transition_sn):
    con = sqlite3.connect(infile)
    
    if not check_sqlite_table(con, "SCORE_MS2"):
        raise click.ClickException("Transition-level scoring for glycoform inference requires prior MS2 or MS1MS2-level scoring. Please run 'pyprophet score --level=ms2' or 'pyprophet score --level=ms1ms2' on this file first.")
    if not check_sqlite_table(con, "FEATURE_TRANSITION"):
        raise click.ClickException("Transition-level feature table not present in file.")

    con.executescript('''
CREATE INDEX IF NOT EXISTS idx_transition_id ON TRANSITION (ID);
CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_precursor_precursor_id ON PRECURSOR (ID);
CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);
CREATE INDEX IF NOT EXISTS idx_feature_feature_id ON FEATURE (ID);
CREATE INDEX IF NOT EXISTS idx_feature_transition_feature_id ON FEATURE_TRANSITION (FEATURE_ID);
CREATE INDEX IF NOT EXISTS idx_feature_transition_transition_id ON FEATURE_TRANSITION (TRANSITION_ID);
''')
    
    table = pd.read_sql_query('''
SELECT TRANSITION.DECOY AS DECOY,
       FEATURE_TRANSITION.*,
       PRECURSOR.CHARGE AS PRECURSOR_CHARGE,
       TRANSITION.PRODUCT_CHARGE AS PRODUCT_CHARGE,
       RUN_ID || '_' || FEATURE_TRANSITION.FEATURE_ID || '_' || PRECURSOR_ID || '_' || TRANSITION_ID AS GROUP_ID
FROM FEATURE_TRANSITION
INNER JOIN
  (SELECT RUN_ID,
          ID,
          PRECURSOR_ID,
          EXP_RT
   FROM FEATURE) AS FEATURE ON FEATURE_TRANSITION.FEATURE_ID = FEATURE.ID
INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID
INNER JOIN
  (SELECT ID,
          CHARGE AS PRODUCT_CHARGE,
          DECOY
   FROM TRANSITION) AS TRANSITION ON FEATURE_TRANSITION.TRANSITION_ID = TRANSITION.ID
WHERE RANK <= %s
  AND PEP <= %s
  AND VAR_ISOTOPE_OVERLAP_SCORE <= %s
  AND VAR_LOG_SN_SCORE > %s
  AND PRECURSOR.DECOY == 0
ORDER BY RUN_ID,
         PRECURSOR.ID,
         FEATURE.EXP_RT,
         TRANSITION.ID;
''' % (max_peakgroup_rank, max_peakgroup_pep, max_transition_isotope_overlap, min_transition_sn), con)

    table.columns = [col.lower() for col in table.columns]
    
    return table


class GlycoformProphetRunner(ProphetRunnerBase):
    def __init__(self, infile, outfile, classifier, 
                 xgb_hyperparams, xgb_params, xgb_params_space, 
                 xeval_fraction, xeval_num_iter, 
                 ss_initial_fdr, ss_iteration_fdr, ss_num_iter, ss_main_score, 
                 group_id, 
                 parametric, pfdr, 
                 pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0, 
                 lfdr_truncate, lfdr_monotone, lfdr_transformation, 
                 lfdr_adj, lfdr_eps, 
                 level, 
                 max_peakgroup_rank, max_peakgroup_pep, 
                 max_transition_isotope_overlap, min_transition_sn,
                 tric_chromprob, threads, test):        
        self.max_peakgroup_rank = max_peakgroup_rank
        self.max_peakgroup_pep = max_peakgroup_pep
        self.max_transition_isotope_overlap = max_transition_isotope_overlap
        self.min_transition_sn = min_transition_sn
        
        super(GlycoformProphetRunner, self).__init__(
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
        self.parametric = parametric
        self.pfdr = pfdr
        self.pi0_lambda = pi0_lambda
        self.pi0_method = pi0_method
        self.pi0_smooth_df = pi0_smooth_df
        self.pi0_smooth_log_pi0 = pi0_smooth_log_pi0        
        self.lfdr_truncate = lfdr_truncate
        self.lfdr_monotone = lfdr_monotone
        self.lfdr_transformation = lfdr_transformation
        self.lfdr_adj = lfdr_adj
        self.lfdr_eps = lfdr_eps
        
        self.tric_chromprob = tric_chromprob
        self.threads = threads
        self.test = test
        
        
    def load_osw(self):
        if self.level == 'ms1':
            table = read_osw_ms1_feature(
                self.infile, 
                max_peakgroup_rank=self.max_peakgroup_rank
            )
        elif self.level == 'transition':
            table = read_osw_transition_feature(
                self.infile, 
                max_peakgroup_rank=self.max_peakgroup_rank,
                max_peakgroup_pep=self.max_peakgroup_pep,
                max_transition_isotope_overlap=self.max_transition_isotope_overlap,
                min_transition_sn=self.min_transition_sn
            )
        else:
            raise click.ClickException("Unspecified data level selected.")
        return table
        
    
    def score(self):
        start_at_glycan = time.time()
        click.echo("Info: scoring glycan part")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            (result, scorer, weights) = self.pyprophet_score()
        
        end_at_glycan = time.time()
        needed = end_at_glycan - start_at_glycan
        seconds = int(needed)
        msecs = int(1000 * (needed - seconds))
        click.echo("Info: glycan part scored: %d seconds and %d msecs" % (seconds, msecs))
    
        result = {
            'scored_table': result.scored_tables,
            'final_statistics': result.final_statistics,
            'summary_statistics': result.summary_statistics
        }
        pi0 = scorer.pi0
        
        self.print_summary(result)
        
        if isinstance(weights, pd.DataFrame):
            weights = weights.assign(part='glycoform')
        else:
            weights = {
               'glycoform': weights
            }
            
        return result, weights, pi0
            
    
    def pyprophet_score(self):
        raise NotImplementedError('pyprophet_score')
        
        
    def extra_writes(self):
        d = super(GlycoformProphetRunner, self).extra_writes()
        d.update({
            'report_path': os.path.join(self.prefix + '_report.pdf')
        })
        return d
    
    
    def save_tsv_results(self, result, extra_writes, pi0):
        super(GlycoformProphetRunner, self) \
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

        df = result['scored_table']
        
        if self.level == "ms1":
            c = con.cursor()
            c.execute('DROP TABLE IF EXISTS SCORE_MS1;')
            con.commit()
            c.fetchall()
            
            table = "SCORE_MS1"
            if 'h_score' in df.columns:
                df = df[['feature_id','d_score','h_score','h0_score','peak_group_rank','p_value','q_value','pep']]
                df.columns = ['FEATURE_ID','SCORE','HSCORE','H0SCORE','RANK','PVALUE','QVALUE','PEP']
            else:
                df = df[['feature_id','d_score','peak_group_rank','p_value','q_value','pep']]
                df.columns = ['FEATURE_ID','SCORE','RANK','PVALUE','QVALUE','PEP']
                
        elif self.level == "transition":
            c = con.cursor()
            c.execute('DROP TABLE IF EXISTS SCORE_TRANSITION;')
            con.commit()
            c.fetchall()
            
            table = "SCORE_TRANSITION"
            df = df[['feature_id','transition_id','d_score','peak_group_rank','p_value','q_value','pep']]
            df.columns = ['FEATURE_ID','TRANSITION_ID','SCORE','RANK','PVALUE','QVALUE','PEP']
        
        df.to_sql(table, con, index=False)        
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
      
    
class GlycoformProphetLearner(GlycoformProphetRunner):
    def pyprophet_score(self):
        table = self.table
        if self.level == 'ms1':
            if 'decoy' in table.columns:
                table = table.drop(columns=['decoy'])
            table = table.rename(columns={'decoy_glycan': 'decoy'})
        
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
            threads=self.threads, test=self.test
        ).learn_and_apply(table)
        
        return (result, scorer, weights)
    
    
    def extra_writes(self):
        d = super(GlycoformProphetLearner, self).extra_writes()
        d.update({
            'trained_weights_path': os.path.join(self.prefix + '_weights.csv'),
            'trained_model_path_ms1': os.path.join(self.prefix + '_ms1_model.bin'),
            'trained_model_path_transition': os.path.join(self.prefix + '_transition_model.bin')
        })
        return d
                
    
class GlycoformProphetWeightApplier(GlycoformProphetRunner):
    def __init__(self, apply_weights, **kwargs):
        super(GlycoformProphetWeightApplier, self) \
            .__init__(**kwargs)
        
        self.persisted_weights = self.load_weights(apply_weights)
    
    
    def pyprophet_score(self):
        if isinstance(self.persisted_weights, pd.DataFrame):
            weights = self.persisted_weights \
                .loc[self.persisted_weights['part'] == 'glycoform'] \
                .drop(columns=['part'])
        else:
            weights = self.persisted_weights['glycoform']
        
        table = self.table
        if self.level == 'ms1':
            if 'decoy' in table.columns:
                table = table.drop(columns=['decoy'])
            table = table.rename(columns={'decoy_glycan': 'decoy'})
        
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
            threads=self.threads, test=self.test
        ).apply_weights(table, weights)
        
        return (result, scorer, weights)    

        
        
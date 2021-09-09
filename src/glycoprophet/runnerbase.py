import click
import os
import pandas as pd
import sqlite3
import pickle

from pyprophet.data_handling import is_sqlite_file


class ProphetRunnerBase:
    def __init__(self, infile, outfile, classifier, 
                 ss_main_score, level):
        self.infile = infile
        self.outfile = outfile
        self.classifier = classifier
        self.ss_main_score = ss_main_score
        self.level = level
        
        if is_sqlite_file(infile):
            self.mode = 'osw'
            self.table = self.load_osw()
            self.table = self.prepare_table(self.table)
        else:
            self.mode = 'tsv'
            self.table = self.load_csv()

        self.prefix = os.path.splitext(outfile)[0]
        
        
    def load_tsv(self):
        table = pd.read_csv(self.infile, "\t")
        return table
    
    def load_osw(self):
        raise NotImplementedError('load_osw')    
        
    def prepare_table(self, table):        
        # Mark main score column
        if self.ss_main_score.lower() in table.columns:
            table = table.rename(index=str, columns={self.ss_main_score.lower(): "main_"+self.ss_main_score.lower()})
        elif self.ss_main_score.lower() == "swath_pretrained":
            # Add a pretrained main score corresponding to the original implementation in OpenSWATH
            # This is optimized for 32-windows SCIEX TripleTOF 5600 data
            table['main_var_pretrained'] = -( -0.19011762 * table['var_library_corr']
                                            +  2.47298914 * table['var_library_rmsd']
                                            +  5.63906731 * table['var_norm_rt_score']
                                            + -0.62640133 * table['var_isotope_correlation_score']
                                            +  0.36006925 * table['var_isotope_overlap_score']
                                            +  0.08814003 * table['var_massdev_score']
                                            +  0.13978311 * table['var_xcorr_coelution']
                                            + -1.16475032 * table['var_xcorr_shape']
                                            + -0.19267813 * table['var_yseries_score']
                                            + -0.61712054 * table['var_log_sn_score'])
        else:
            raise click.ClickException("Main score column not present in data.")
    
        # Enable transition count & precursor / product charge scores for XGBoost-based classifier
        if self.classifier == 'XGBoost':
            click.echo("Info: Enable number of transitions & precursor / product charge scores for XGBoost-based classifier")
            table = table.rename(index=str, columns={'precursor_charge': 'var_precursor_charge', 'product_charge': 'var_product_charge', 'transition_count': 'var_transition_count'})
        
        return table
    
    
    def score(self):
        raise NotImplementedError('score')
        
    def print_summary(self, result):
        summary = result.get('summary_statistics', None)  
        if summary is not None:
            click.echo("=" * 80)
            click.echo(summary)
            click.echo("=" * 80)

    
    def extra_writes(self):
        return {
            'output_path': os.path.join(self.prefix + '_scored.tsv'),
            'summ_stat_path': os.path.join(self.prefix + '_summary_stat.csv'),
            'full_stat_path': os.path.join(self.prefix + '_full_stat.csv')
        }
        
    
    def save_osw_results(self, result, extra_writes, pi0):
        raise NotImplementedError('save_osw_results')
        
    
    def save_tsv_results(self, result, extra_writes, pi0):
        summ_stat_path = extra_writes.get("summ_stat_path", None)
        if summ_stat_path is not None:
            result['summary_statistics'].to_csv(summ_stat_path, sep=",", index=False)
        click.echo("Info: %s written." % summ_stat_path)

        full_stat_path = extra_writes.get("full_stat_path", None)
        if full_stat_path is not None:
            result['final_statistics'].to_csv(full_stat_path, sep=",", index=False)
            click.echo("Info: %s written." % full_stat_path)

        output_path = extra_writes.get("output_path", None)
        if output_path is not None:
            result['scored_table'].to_csv(output_path, sep="\t", index=False)
            click.echo("Info: %s written." % output_path)

    
    def save_osw_weights(self, weights):
        if self.classifier == "LDA":
            weights['level'] = self.level
            con = sqlite3.connect(self.outfile)

            c = con.cursor()
            c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="GLYCOPEPTIDEPROPHET_WEIGHTS";')
            if c.fetchone()[0] == 1:
                c.execute('DELETE FROM GLYCOPEPTIDEPROPHET_WEIGHTS WHERE LEVEL =="%s"' % self.level)
            c.close()

            weights.to_sql("GLYCOPEPTIDEPROPHET_WEIGHTS", con, index=False, if_exists='append')

        elif self.classifier == "XGBoost":
            con = sqlite3.connect(self.outfile)

            c = con.cursor()
            c.execute('SELECT count(name) FROM sqlite_master WHERE type="table" AND name="GLYCOPEPTIDEPROPHET_XGB";')
            if c.fetchone()[0] == 1:
                c.execute('DELETE FROM GLYCOPEPTIDEPROPHET_XGB WHERE LEVEL =="%s"' % self.level)
            else:
                c.execute('CREATE TABLE GLYCOPEPTIDEPROPHET_XGB (level TEXT, xgb BLOB)')

            c.execute('INSERT INTO GLYCOPEPTIDEPROPHET_XGB VALUES(?, ?)', [self.level, pickle.dumps(weights)])
            con.commit()
            c.close()
            
        
    def save_bin_weights(self, weights, extra_writes):
        trained_weights_path = extra_writes.get("trained_model_path_" + self.level)
        if trained_weights_path is not None:
            with open(trained_weights_path, 'wb') as file:
                self.persisted_weights = pickle.dump(weights, file)
            click.echo("Info: %s written." % trained_weights_path)
    
    
    def save_tsv_weights(self, weights, extra_writes):
        weights['level'] = self.level
        trained_weights_path = extra_writes.get("trained_weights_path")
        if trained_weights_path is not None:
            weights.to_csv(trained_weights_path, sep=",", index=False)
            click.echo("Info: %s written." % trained_weights_path)
            
    
    def load_tsv_weights(self, apply_weights):
        try:
            persisted_weights = pd.read_csv(apply_weights, sep=",")
            if self.level != persisted_weights['level'].unique()[0]:
                raise click.ClickException("Weights file has wrong level.")
        except Exception:
            import traceback
            traceback.print_exc()
            raise
        return persisted_weights
    
    
    def load_bin_weights(self, apply_weights):
        with open(apply_weights, 'rb') as file:
            persisted_weights = pickle.load(file)            
        return persisted_weights
            
    
    def load_osw_weights(self, apply_weights):
        if self.classifier == "LDA":
            try:
                con = sqlite3.connect(apply_weights)
                data = pd.read_sql_query("SELECT * FROM GLYCOPEPTIDEPROPHET_WEIGHTS WHERE LEVEL=='%s'" % self.level, con)
                data.columns = [col.lower() for col in data.columns]
                con.close()
                persisted_weights = data
                if self.level != persisted_weights['level'].unique()[0]:
                    raise click.ClickException("Weights file has wrong level.")
            except Exception:
                import traceback
                traceback.print_exc()
                raise
        elif self.classifier == "XGBoost":
            try:
                con = sqlite3.connect(apply_weights)
                data = con.execute("SELECT xgb FROM PYPROPHET_XGB WHERE LEVEL=='%s'" % self.level).fetchone()
                con.close()
                persisted_weights = pickle.loads(data[0])
            except Exception:
                import traceback
                traceback.print_exc()
                raise
        return persisted_weights
    
    
    def load_weights(self, apply_weights):
        if not os.path.exists(apply_weights):
            raise click.ClickException("Weights file %s does not exist." % apply_weights)
        
        if self.mode == "tsv":
            if self.classifier == "LDA":
                return self.load_tsv_weights(apply_weights)
            elif self.classifier == "XGBoost":
                return self.load_bin_weights(apply_weights)
        elif self.mode == "osw":
            return self.load_osw_weights(apply_weights)
        
    
    def run(self):
        extra_writes = dict(self.extra_writes())
          
        result, weights, pi0 = self.score()
        
        if self.mode == 'tsv':
            self.save_tsv_results(result, extra_writes, pi0)
            if self.classifier == 'LDA':
                self.save_tsv_weights(weights, extra_writes)
            elif self.classifier == 'XGBoost':
                self.save_bin_weights(weights, extra_writes)

        elif self.mode == 'osw':
            self.save_osw_results(result, extra_writes, pi0)
            self.save_osw_weights(weights)
            
            
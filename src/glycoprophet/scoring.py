import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

class LDAGlycoPeptideScorer():
    def __init__(self, group_id='group_id', score='d_score_combined'):
        self.group_id = group_id
        self.score = score
    
    def rerank_table(self, table, ascending=False, **kwargs):
        def rerank_peak_groups(table):        
            table = table.sort_values(by=self.score, ascending=ascending, **kwargs)
            table['peak_group_rank'] = list(range(1, len(table) + 1))
            return table
        
        table = table.groupby(self.group_id) \
            .apply(rerank_peak_groups)
        table.reset_index(drop=True, inplace=True)
        return table

    
    def learn_and_apply(self, table):
        if 'peak_group_rank' not in table.columns:
            if self.score not in table.columns:
                table[self.score] = \
                    0.5 * table['d_score_peptide'] + \
                    0.5 * table['d_score_glycan']
            table = self.rerank_table(table)
        
        table_learning = table.loc[table['peak_group_rank'] == 1]
        
        classifier = LinearDiscriminantAnalysis()
        classifier.fit(
            table_learning[['d_score_peptide', 'd_score_glycan']], 
            table_learning['decoy_glycan'] | table_learning['decoy_peptide']
        )
        coef = classifier.scalings_.flatten() / \
            np.sum(classifier.scalings_.flatten())
        
        score = \
            table['d_score_peptide'] * coef[0] + \
            table['d_score_glycan'] * coef[1]
        
        table[self.score] = score
        table = self.rerank_table(table)
        weights = pd.DataFrame.from_dict({
            'score': ['d_score_peptide', 'd_score_glycan'],
            'weight': coef
        })
        return {
            'scored_table': table,
            'weights': weights
        }
        
    
    def apply_weights(self, table, loaded_weights):
        coef = np.array([
            loaded_weights \
                .loc[loaded_weights['score'] == 'd_score_' + x, 'weight'] \
                .values[0]
            for x in ['peptide', 'glycan']
        ])
        
        score = \
            table['d_score_peptide'] * coef[0] + \
            table['d_score_glycan'] * coef[1]
        
        table[self.score] = score
        table = self.rerank_table(table)
        return {
            'scored_table': table,
            'weights': loaded_weights
        }
        
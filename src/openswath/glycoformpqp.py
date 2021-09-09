import pandas as pd
import numpy as np

from .glycopqp import GlycoPeptideQueryParameter
from pepmass.glycomass import GlycanNode

class GlycoformUisQueryParameter(GlycoPeptideQueryParameter):
    def prepare_tables(self, data):
        tables = super(GlycoformUisQueryParameter, self) \
            .prepare_tables(data)
            
        data_iden = data \
            .loc[data['IdentifyingTransition'] != 0]
                
        glycoform_glycopeptide = data_iden \
            [['ModifiedPeptideSequence', 'GlycanSite', 'Glycoform']] \
            .drop_duplicates()
        glycoform_glycopeptide['Glycoform'] = \
            glycoform_glycopeptide['Glycoform'].str.split('|') 
        glycoform_glycopeptide['num_glycoforms'] = \
            glycoform_glycopeptide['Glycoform'].map(len)
        glycoform_glycopeptide = pd.DataFrame({
            'MODIFIED_SEQUENCE': np.repeat(
                glycoform_glycopeptide['ModifiedPeptideSequence'].values,
                glycoform_glycopeptide['num_glycoforms'].values
            ), 
            'GLYCAN_SITE': np.repeat(
                glycoform_glycopeptide['GlycanSite'].values,
                glycoform_glycopeptide['num_glycoforms'].values
            ), 
            'GLYCAN_NAME' : np.concatenate(
                glycoform_glycopeptide['Glycoform'].values
            )
        }) 
        glycoform_glycopeptide.drop_duplicates(inplace=True) 
              
        glycoform_glycopeptide = glycoform_glycopeptide.merge(
            tables['peptide'][['ID', 'MODIFIED_SEQUENCE']] \
                .rename(columns={'ID': 'PEPTIDE_ID'}, copy=False),
            how='inner',
            copy=False
        )
        glycoform_glycopeptide = glycoform_glycopeptide.merge(
            tables['glycan'][['ID', 'GLYCAN_NAME']] \
                .rename(columns={'ID': 'GLYCAN_ID'}, copy=False),
            how='left',
            copy=False            
        )
        glycoform_glycopeptide = glycoform_glycopeptide.merge(
            tables['glycopeptide'] \
                [['ID', 'GLYCAN_SITE', 'PEPTIDE_ID', 'GLYCAN_ID']] \
                .rename(columns={'ID': 'GLYCOPEPTIDE_ID'}, copy=False),
            how='left',
            copy=False            
        )
    
        new_glycan = glycoform_glycopeptide['GLYCAN_ID'].isna()   
        glycan = glycoform_glycopeptide \
            .loc[new_glycan, ['GLYCAN_NAME']] \
            .drop_duplicates()
        max_glycan_id = tables['glycan']['ID'].max()
        if not max_glycan_id >= 0:
            max_glycan_id = 0
        glycan.insert(0, 'ID', np.array(
            range(max_glycan_id + 1, \
                       max_glycan_id + 1 + len(glycan)),
            dtype=np.int64
        ))
        glycan['GLYCAN_STRUCT'] = glycan['GLYCAN_NAME'] \
            .str.replace('DECOY_', '')
        glycan['GLYCAN_COMPOSITION'] = glycan['GLYCAN_STRUCT'] \
            .map(lambda x: GlycanNode.from_str(x).composition_str())
        glycan['DECOY'] = glycan['GLYCAN_NAME'] \
            .str.startswith('DECOY').astype(np.int64)
        glycoform_glycopeptide.loc[new_glycan, 'GLYCAN_ID'] = \
            glycoform_glycopeptide.loc[new_glycan, ['GLYCAN_NAME']] \
                .merge(
                    glycan[['ID', 'GLYCAN_NAME']],
                    on=['GLYCAN_NAME'],
                    how='left',
                    copy=False
                )[['ID']].values        
        glycoform_glycopeptide['GLYCAN_ID'] = \
            glycoform_glycopeptide['GLYCAN_ID'].astype(int)
                
        tables['glycan'] = pd.concat(
            [
                tables['glycan'], 
                glycan[[
                    'ID', 'GLYCAN_NAME', 
                    'GLYCAN_STRUCT', 'GLYCAN_COMPOSITION', 
                    'DECOY'
                ]]
            ], 
            axis=0, 
            ignore_index=True
        )
        
        glycoform_glycopeptide['GLYCOPEPTIDE_NAME'] = \
            glycoform_glycopeptide.agg(lambda x: \
                x['MODIFIED_SEQUENCE'] + \
                '_' + str(x['GLYCAN_SITE']) + \
                ',' + str(x['GLYCAN_NAME']), \
                axis=1
            )
        
        new_glycopeptide = glycoform_glycopeptide['GLYCOPEPTIDE_ID'].isna()   
        glycopeptide = glycoform_glycopeptide \
            .loc[new_glycopeptide][[
                'GLYCOPEPTIDE_NAME', 
                'GLYCAN_SITE', 
                'PEPTIDE_ID',
                'GLYCAN_ID'
            ]] \
            .drop_duplicates()
        max_glycopeptide_id = tables['glycopeptide']['ID'].max()
        if not max_glycopeptide_id >= 0:
            max_glycopeptide_id = 0
        glycopeptide.insert(0, 'ID', np.array(
            range(max_glycopeptide_id + 1, \
                       max_glycopeptide_id + 1 + len(glycopeptide)),
            dtype=np.int64
        ))
        glycopeptide.insert(
            3, 'DECOY_PEPTIDE', 
            glycopeptide[['PEPTIDE_ID']].merge(
                tables['peptide'][['ID', 'DECOY']] \
                    .rename(columns={'ID': 'PEPTIDE_ID'}, copy=False),
                how='left',
                copy=False
            )['DECOY'].values
        )
        glycopeptide.insert(
            4, 'DECOY_GLYCAN', 
            glycopeptide[['GLYCAN_ID']].merge(
                tables['glycan'][['ID', 'DECOY']]
                    .rename(columns={'ID': 'GLYCAN_ID'}, copy=False),
                how='left',
                copy=False
            )['DECOY'].values
        )
        glycopeptide['GLYCOSITE_ID'] = \
            glycopeptide[['GLYCAN_SITE', 'PEPTIDE_ID']].merge(
                tables['glycopeptide'] \
                    [['GLYCAN_SITE', 'PEPTIDE_ID', 'GLYCOSITE_ID']] \
                    .drop_duplicates(subset=['GLYCAN_SITE', 'PEPTIDE_ID']),
                how='left',
                copy=False
            )['GLYCOSITE_ID'].values        
        glycoform_glycopeptide.loc[new_glycopeptide, 'GLYCOPEPTIDE_ID'] = \
            glycoform_glycopeptide \
                .loc[new_glycopeptide, ['GLYCOPEPTIDE_NAME']] \
                .merge(
                    glycopeptide[['ID', 'GLYCOPEPTIDE_NAME']],
                    on=['GLYCOPEPTIDE_NAME'],
                    how='left',
                    copy=False
                )['ID'].values
        glycoform_glycopeptide['GLYCOPEPTIDE_ID'] = \
            glycoform_glycopeptide['GLYCOPEPTIDE_ID'].astype(np.int64)
        
        tables['glycopeptide'] = pd.concat(
            [
                tables['glycopeptide'], 
                glycopeptide[[
                    'ID', 'GLYCOPEPTIDE_NAME', 'GLYCAN_SITE',
                    'DECOY_PEPTIDE', 'DECOY_GLYCAN',
                    'PEPTIDE_ID', 'GLYCAN_ID', 'GLYCOSITE_ID'
                ]]        
            ], 
            axis=0, 
            ignore_index=True
        ) 
        
        transition_glycopeptide = data_iden[[
            'TransitionId', 
            'ModifiedPeptideSequence', 
            'GlycanSite',
            'Glycoform'
        ]].rename(columns={
            'TransitionId': 'TRAML_ID', 
            'ModifiedPeptideSequence': 'MODIFIED_SEQUENCE', 
            'GlycanSite': 'GLYCAN_SITE'
        }, copy=False)
        transition_glycopeptide = transition_glycopeptide.merge(
            tables['transition'][['ID', 'TRAML_ID']],
            how='inner',
            copy=False
        )        
            
        
        def match_transition_glycopeptide(transition_glycopeptide):
            glycoform = \
                transition_glycopeptide['Glycoform'].str.split('|')             
            num_glycoforms = \
                glycoform.map(len)
            transition_glycopeptide = pd.DataFrame({
                'TRANSITION_ID': np.repeat(
                    transition_glycopeptide['ID'].values,
                    num_glycoforms.values
                ), 
                'MODIFIED_SEQUENCE': np.repeat(
                    transition_glycopeptide['MODIFIED_SEQUENCE'].values,
                    num_glycoforms.values
                ),
                'GLYCAN_SITE': np.repeat(
                    transition_glycopeptide['GLYCAN_SITE'].values,
                    num_glycoforms.values
                ), 
                'GLYCAN_NAME' : np.concatenate(
                    glycoform.values
                )
            }) 
            
            transition_glycopeptide = transition_glycopeptide.merge(
                glycoform_glycopeptide[[
                    'MODIFIED_SEQUENCE', 'GLYCAN_SITE',
                    'GLYCAN_NAME', 'GLYCOPEPTIDE_ID'
                ]],
                on=['MODIFIED_SEQUENCE', 'GLYCAN_SITE', 'GLYCAN_NAME'],
                how='inner',
                copy=False
            )[['TRANSITION_ID', 'GLYCOPEPTIDE_ID']]
            transition_glycopeptide.sort_values('TRANSITION_ID', inplace=True)
            
            return transition_glycopeptide
        
        subset_size = 100000
        
        transition_glycopeptide = pd.concat((
            match_transition_glycopeptide(
                transition_glycopeptide \
                    .iloc[i:(i + subset_size)]
            )
            for i in range(0, len(transition_glycopeptide), subset_size)            
        ), axis=0, ignore_index=True)
                               
        tables['transition_glycopeptide_mapping'] = transition_glycopeptide[[
            'TRANSITION_ID', 'GLYCOPEPTIDE_ID'
        ]]       
        
        return tables
    
        
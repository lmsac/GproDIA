from .pqp import PeptideQueryParameter
from . import pqp_sqlconst as sqlconst
from . import glycopqp_sqlconst as gsqlconst
import pandas as pd

class GlycoPeptideQueryParameter(PeptideQueryParameter):
    def __init__(self, path):
        super(GlycoPeptideQueryParameter, self).__init__(path)
        
    def create_table(self, skip_if_exist=False): 
        super(GlycoPeptideQueryParameter, self) \
            .create_table(skip_if_exist=skip_if_exist)
        
        SQL_CREATE_TABLE = [
            getattr(gsqlconst, x)
            for x in dir(gsqlconst)
            if x.startswith('CREATE_TABLE_')
        ]
        
        if skip_if_exist:
            SQL_CREATE_TABLE = [
                x.replace('CREATE TABLE ', 'CREATE TABLE IF NOT EXISTS ')
                for x in SQL_CREATE_TABLE
            ]
            
        try:
            cursor = self.connection.cursor()            
            for sql in SQL_CREATE_TABLE:
                cursor.execute(sql)        
            
            self.connection.commit()        
        except Exception:
            self.connection.rollback()
            raise
        
    
    def get_glycan(self, glycan_id=None, 
                   glycan_name=None,
                   glycan_struct=None,
                   glycan_composition=None):
        sql = gsqlconst.SELECT_GLYCAN
        
        where = []
        params = {}
        if glycan_id is not None:
            where.append('ID = :glycan_id')
            params['glycan_id'] = glycan_id
        if glycan_name is not None:
            where.append('GLYCAN_NAME = :glycan_name')
            params['glycan_name'] = glycan_name        
        if glycan_struct is not None:
            where.append('GLYCAN_STRUCT = :glycan_struct')
            params['glycan_struct'] = glycan_struct
        if glycan_composition is not None:
            where.append('GLYCAN_COMPOSITION = :glycan_composition')
            params['glycan_composition'] = glycan_composition
            
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
    
        
    def get_glycosite(self, glycosite_id=None, 
                      glycosite_name=None,
                      protein_id=None, include_protein=None):   
        if protein_id is not None:
            include_protein = True    
            
        sql = gsqlconst.SELECT_GLYCOSITE          
        if include_protein:
            sql += gsqlconst.JOIN_GLYCOSITE_PROTEIN
            
        where = []
        params = {}
        if glycosite_id is not None:
            where.append('ID = :glycosite_id')
            params['glycosite_id'] = glycosite_id
        if glycosite_name is not None:
            where.append('GLYCOSITE_NAME = :glycosite_name')
            params['glycosite_name'] = glycosite_name        
        if protein_id is not None:
            where.append('PROTEIN_ID = :protein_id')
            params['protein_id'] = protein_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
        
            
    def get_glycopeptide(self, glycopeptide_id=None, 
                         glycopeptide_name=None, glycan_site=None,
                         peptide_id=None, include_peptide=None,
                         glycosite_id=None, include_glycosite=None,
                         glycan_id=None, include_glycan=None):        
        if peptide_id is not None:
            include_peptide = True
        if glycosite_id is not None:
            include_glycosite = True
        if glycan_id is not None:
            include_glycan = True
        
        sql = gsqlconst.SELECT_GLYCOPEPTIDE              
        if include_peptide:
            sql += gsqlconst.JOIN_GLYCOPEPTIDE_PEPTIDE    
        if include_glycosite:
            sql += gsqlconst.JOIN_GLYCOPEPTIDE_GLYCOSITE
        if include_glycan:
            sql += gsqlconst.JOIN_GLYCOPEPTIDE_GLYCAN
        
        where = []
        params = {}
        if glycopeptide_id is not None:
            where.append('ID = :glycopeptide_id')
            params['glycopeptide_id'] = glycopeptide_id
        if glycopeptide_name is not None:
            where.append('GLYCOPEPTIDE_NAME = :glycopeptide_name')
            params['glycopeptide_name'] = glycopeptide_name        
        if glycan_site is not None:
            where.append('GLYCAN_SITE = :glycan_site')
            params['glycan_site'] = glycan_site
        if peptide_id is not None:
            where.append('PEPTIDE_ID = :peptide_id')
            params['peptide_id'] = peptide_id
        if glycosite_id is not None:
            where.append('GLYCOSITE_ID = :glycosite_id')
            params['glycosite_id'] = glycosite_id
        if glycan_id is not None:
            where.append('GLYCAN_ID = :glycan_id')
            params['glycan_id'] = glycan_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
    
    
    def get_precursor(self, precursor_id=None,
                      precursor_name=None,
                      glycopeptide_id=None, include_glycopeptide=None,
                      peptide_id=None, include_peptide=None,
                      compound_id=None, include_compound=None):        
        if glycopeptide_id is not None:
            include_glycopeptide = True
        if peptide_id is not None:
            include_peptide = True
        if compound_id is not None:
            include_compound = True
        
        sql = sqlconst.SELECT_PRECURSOR
        if include_glycopeptide:
            sql += gsqlconst.JOIN_PRECURSOR_GLYCOPEPTIDE
        if include_peptide:
            sql += sqlconst.JOIN_PRECURSOR_PEPTIDE
        if include_compound:
            sql += sqlconst.JOIN_PRECURSOR_COMPOUND
            
        where = []
        params = {}
        if precursor_id is not None:
            where.append('ID = :precursor_id')
            params['precursor_id'] = precursor_id
        if precursor_name is not None:
            where.append('TRAML_ID = :precursor_name')
            params['precursor_name'] = precursor_name
        if glycopeptide_id is not None:
            where.append('GLYCOPEPTIDE_ID = :glycopeptide_id')
            params['glycopeptide_id'] = glycopeptide_id
        if peptide_id is not None:
            where.append('PEPTIDE_ID = :peptide_id')
            params['peptide_id'] = peptide_id
        if compound_id is not None:
            where.append('COMPOUND_ID = :compound_id')
            params['compound_id'] = compound_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
    
    
    def add_glycan(self, glycan): 
        self.add_table(glycan, 'GLYCAN')
        
    
    def add_glycosite(self, glycosite): 
        if 'PROTEIN_ID' in glycosite.columns:
            glycosite_protein_mapping = \
                glycosite[['ID', 'PROTEIN_ID']] \
                .rename(columns={'ID': 'GLYCOSITE_ID'})   
            glycosite = glycosite.drop(columns=['PROTEIN_ID'])
        else:
            glycosite_protein_mapping
            
        self.add_table(glycosite, 'GLYCOSITE')
        
        if glycosite_protein_mapping is not None:
            self.add_table(
                glycosite_protein_mapping,
                'GLYCOSITE_PROTEIN_MAPPING'
            )
    
    
    def add_glycopeptide(self, glycopeptide):
        if 'GLYCOSITE_ID' in glycopeptide.columns:
            glycopeptide_glycosite_mapping = \
                glycopeptide[['ID', 'GLYCOSITE_ID']] \
                .rename(columns={'ID': 'GLYCOPEPTIDE_ID'})   
            glycopeptide = glycopeptide.drop(columns=['GLYCOSITE_ID'])
        else:
            glycopeptide_glycosite_mapping = None
        if 'GLYCAN_ID' in glycopeptide.columns:
            glycopeptide_glycan_mapping = \
                glycopeptide[['ID', 'GLYCAN_ID']] \
                .rename(columns={'ID': 'GLYCOPEPTIDE_ID'})   
            glycopeptide = glycopeptide.drop(columns=['GLYCAN_ID'])
        else:
            glycopeptide_glycan_mapping = None
        if 'PEPTIDE_ID' in glycopeptide.columns:
            glycopeptide_peptide_mapping = \
                glycopeptide[['ID', 'PEPTIDE_ID']] \
                .rename(columns={'ID': 'GLYCOPEPTIDE_ID'})   
            glycopeptide = glycopeptide.drop(columns=['PEPTIDE_ID'])
        else:
            glycopeptide_peptide_mapping = None
            
        self.add_table(glycopeptide, 'GLYCOPEPTIDE')
        
        if glycopeptide_glycosite_mapping is not None:
            self.add_table(
                glycopeptide_glycosite_mapping,
                'GLYCOPEPTIDE_GLYCOSITE_MAPPING'
            )
        if glycopeptide_glycan_mapping is not None:
            self.add_table(
                glycopeptide_glycan_mapping,                           
                'GLYCOPEPTIDE_GLYCAN_MAPPING'
            )
        if glycopeptide_peptide_mapping is not None:
            self.add_table(
                glycopeptide_peptide_mapping,
                'GLYCOPEPTIDE_PEPTIDE_MAPPING'
            )
        
        
        
    def add_precursor(self, precursor):                
        if 'GLYCOPEPTIDE_ID' in precursor.columns:
            precursor_glycopeptide_mapping = precursor[['ID', 'GLYCOPEPTIDE_ID']] \
                .rename(columns={'ID': 'PRECURSOR_ID'})   
            precursor = precursor.drop(columns=['GLYCOPEPTIDE_ID'])
        else:
            precursor_glycopeptide_mapping = None            
          
        super(GlycoPeptideQueryParameter, self) \
            .add_precursor(precursor)
        
        if precursor_glycopeptide_mapping is not None:
            self.add_table(
                precursor_glycopeptide_mapping, 
                'PRECURSOR_GLYCOPEPTIDE_MAPPING'
            )
            
    
    def prepare_tables(self, data):
        tables = super(GlycoPeptideQueryParameter, self) \
            .prepare_tables(data)
        
        cursor = self.connection.cursor()
               
        glycan = data[['GlycanId', 'GlycanStruct', 'GlycanComposition', 
                       'DecoyGlycan']] \
            .drop_duplicates(['GlycanId'])
            
        cursor.execute(gsqlconst.SELECT_MAX_ID_GLYCAN)
        glycan_max_id = cursor.fetchone()[0]
        
        if glycan_max_id is None:
            glycan['Exist'] = False
            glycan['ID'] = list(range(0, len(glycan)))
            glycan_max_id = len(glycan)
            
        else:
            glycan_id = [
                next((x for i, x in self \
                      .get_glycan(glycan_struct=x)['ID'] \
                      .iteritems()), None)
                for x in glycan['GlycanId'].values
            ]
            glycan['Exist'] = list(map(lambda x: x is not None, 
                glycan_id))
            for i, x in enumerate(glycan_id):
                if x is None:
                    glycan_max_id += 1
                    glycan_id[i] = glycan_max_id
                    
            glycan['ID'] = glycan_id
            
        glycan.set_index('GlycanId', drop=False, inplace=True)
            
        
        glycosite = data[[
            'ProteinId',
            'GlycoSiteId',
            'ProteinGlycoSite',
            'DecoyPeptide'
        ]].drop_duplicates(['GlycoSiteId'])
        
        cursor.execute(gsqlconst.SELECT_MAX_ID_GLYCOSITE)
        glycosite_max_id = cursor.fetchone()[0]
        
        if glycosite_max_id is None:
            glycosite['Exist'] = False
            glycosite['ID'] = list(range(0, len(glycosite)))
            glycosite_max_id = len(glycosite)
            
        else:
            glycosite_id = [
                next((x for i, x in self \
                      .get_glycosite(glycosite_name=x)['ID'] \
                      .iteritems()), None)
                for x in zip(glycosite['GlycoSiteId'].values)
            ]
            glycosite['Exist'] = list(map(lambda x: x is not None, 
                glycosite_id))
            for i, x in enumerate(glycosite_id):
                if x is None:
                    glycosite_max_id += 1
                    glycosite_id[i] = glycosite_max_id
                    
            glycosite['ID'] = glycosite_id
            
        glycosite.set_index('GlycoSiteId', drop=False, inplace=True)
                 
                
        glycopeptide = data[[
            'ModifiedPeptideSequence', 'GlycanId', 'GlycoSiteId',
            'GlycoPeptideId', 'GlycanSite',
            'DecoyPeptide', 'DecoyGlycan'            
        ]].drop_duplicates(['GlycoPeptideId'])
        
        cursor.execute(gsqlconst.SELECT_MAX_ID_GLYCOPEPTIDE)
        glycopeptide_max_id = cursor.fetchone()[0]
        
        if glycopeptide_max_id is None:
            glycopeptide['Exist'] = False
            glycopeptide['ID'] = list(range(0, len(glycopeptide)))
            glycopeptide_max_id = len(glycopeptide)
        
        else:
            glycopeptide_id = [
                next((x for i, x in self \
                      .get_glycopeptide(glycopeptide_name=x)['ID'] \
                      .iteritems()), None)
                for x in glycopeptide['GlycoPeptideId'].values
            ]
            glycopeptide['Exist'] = list(map(lambda x: x is not None, 
                glycopeptide_id))
            for i, x in enumerate(glycopeptide_id):
                if x is None:
                    glycopeptide_max_id += 1
                    glycopeptide_id[i] = glycopeptide_max_id
                    
            glycopeptide['ID'] = glycopeptide_id
        
        glycopeptide.set_index('GlycoPeptideId', drop=False, inplace=True)
              
        
        precursor = tables['precursor']                  
        precursor = pd.merge(
            precursor,
            data[['TransitionGroupId', 'GlycoPeptideId']].drop_duplicates(),
            how='left', left_on='TRAML_ID', right_on='TransitionGroupId'
        )                     
        precursor['GLYCOPEPTIDE_ID'] = [
            glycopeptide.loc[x['GlycoPeptideId'], 'ID']
            for i, x in precursor.iterrows()
        ]
        precursor.drop(
            columns=['TransitionGroupId', 'GlycoPeptideId'], 
            inplace=True
        )
        
        
        peptide = tables['peptide']
        protein = tables['protein']
                
        glycopeptide = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['GlycoPeptideId'], 
                 x['GlycanSite'],
                 x['DecoyPeptide'], x['DecoyGlycan'],
                 peptide['ID'].iloc[ \
                    peptide['MODIFIED_SEQUENCE'] \
                        .eq(x['ModifiedPeptideSequence']) \
                        .idxmax()],
                 glycan.loc[x['GlycanId'], 'ID'],
                 glycosite.loc[x['GlycoSiteId'], 'ID'])
                for i, x in glycopeptide.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'GLYCOPEPTIDE_NAME',
                'GLYCAN_SITE',
                'DECOY_PEPTIDE',
                'DECOY_GLYCAN',
                'PEPTIDE_ID',
                'GLYCAN_ID',
                'GLYCOSITE_ID'
            ]   
        )
                
        glycosite = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['GlycoSiteId'], 
                 x['ProteinGlycoSite'], x['DecoyPeptide'],
                 protein['ID'].iloc[ \
                    protein['PROTEIN_ACCESSION'] \
                        .eq(x['ProteinId']) \
                        .idxmax()])
                for i, x in glycosite.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'GLYCOSITE_NAME',
                'PROTEIN_GLYCOSITE',
                'DECOY',
                'PROTEIN_ID'                
            ]   
        )
        
        glycan = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['GlycanId'], 
                 x['GlycanStruct'], x['GlycanComposition'],
                 x['DecoyGlycan'])
                for i, x in glycan.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'GLYCAN_NAME',
                'GLYCAN_STRUCT',
                'GLYCAN_COMPOSITION',
                'DECOY'
            ]   
        )
                    
        tables['precursor'] = precursor
        tables['glycopeptide'] = glycopeptide
        tables['glycosite'] = glycosite
        tables['glycan'] = glycan
        
        return tables     
        
            
import sqlite3
import pandas as pd

from . import pqp_sqlconst as sqlconst

class PeptideQueryParameter:
    def __init__(self, path):
        self.connection = sqlite3.connect(path)
        
    
    def close(self):
        if self.connection is not None:
            self.connection.close()
            
    
    def create_table(self, skip_if_exist=False):        
        SQL_CREATE_TABLE = [
            getattr(sqlconst, x)
            for x in dir(sqlconst)
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
        
            cursor.execute(sqlconst.INSERT_VERSION, (3,))
            # cursor.execute(sqlconst.INSERT_GENE, (0, 'NA', 1))
            self.connection.commit()        
        except Exception:
            self.connection.rollback()
            raise
            
            
    def get_gene(self, gene_id=None, gene_name=None):
        sql = sqlconst.SELECT_GENE
        
        where = []
        params = {}
        if gene_id is not None:
            where.append('ID = :gene_id')
            params['gene_id'] = gene_id
        if gene_name is not None:
            where.append('GENE_NAME = :gene_name')
            params['gene_name'] = gene_name
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)


    def get_protein(self, protein_id=None, accession=None):        
        sql = sqlconst.SELECT_PROTEIN
        
        where = []
        params = {}
        if protein_id is not None:
            where.append('ID = :protein_id')
            params['protein_id'] = protein_id
        if accession is not None:
            where.append('PROTEIN_ACCESSION = :accession')
            params['accession'] = accession
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
        
            
    def get_peptide(self, peptide_id=None, 
                    modified_sequence=None, stripped_sequence=None,
                    protein_id=None, include_protein=None,
                    gene_id=None, include_gene=None):
        if protein_id is not None:
            include_protein = True
        if gene_id is not None:
            include_gene = True
        
        sql = sqlconst.SELECT_PEPTIDE
        if include_protein:
            sql += sqlconst.JOIN_PEPTIDE_PROTEIN
        if include_gene:
            sql += sqlconst.JOIN_PEPTIDE_GENE            
            
        where = []
        params = {}
        if peptide_id is not None:
            where.append('ID = :peptide_id')
            params['peptide_id'] = peptide_id
        if modified_sequence is not None:
            where.append('MODIFIED_SEQUENCE = :modified_sequence')
            params['modified_sequence'] = modified_sequence
        if stripped_sequence is not None:
            where.append('UNMODIFIED_SEQUENCE = :stripped_sequence')
            params['stripped_sequence'] = stripped_sequence
        if protein_id is not None:
            where.append('PROTEIN_ID = :protein_id')
            params['protein_id'] = protein_id
        if gene_id is not None:
            where.append('GENE_ID = :gene_id')
            params['gene_id'] = gene_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
        
                    
    def get_precursor(self, precursor_id=None,
                      precursor_name=None,
                      peptide_id=None, include_peptide=None,
                      compound_id=None, include_compound=None):        
        if peptide_id is not None:
            include_peptide = True
        if compound_id is not None:
            include_compound = True
        
        sql = sqlconst.SELECT_PRECURSOR
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
        if peptide_id is not None:
            where.append('PEPTIDE_ID = :peptide_id')
            params['peptide_id'] = peptide_id
        if compound_id is not None:
            where.append('COMPOUND_ID = :compound_id')
            params['compound_id'] = compound_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
            
    
    def get_transition(self, transition_id=None,
                       transition_name=None,
                       precursor_id=None, include_precursor=None):        
        if precursor_id is not None:
            include_precursor = True
        
        sql = sqlconst.SELECT_TRANSITION
        if include_precursor:
            sql += sqlconst.JOIN_TRANSITION_PRECURSOR
            
        where = []
        params = {}
        if transition_id is not None:
            where.append('ID = :transition_id')
            params['transition_id'] = transition_id
        if transition_name is not None:
            where.append('TRAML_ID = :transition_name')
            params['transition_name'] = transition_name
        if precursor_id is not None:
            where.append('PRECURSOR_ID = :precursor_id')
            params['precursor_id'] = precursor_id
        
        if len(where) > 0:
            sql += ' WHERE ' + ' AND '.join(where)
        return pd.read_sql_query(sql, self.connection, params=params)
    
    
    def add_table(self, table, table_name):
        table.to_sql(
            table_name, self.connection,
            if_exists='append', index=False
        )
    
    
    def add_gene(self, gene):
        self.add_table(gene, 'GENE')
        
        
    def add_protein(self, protein):
        self.add_table(protein, 'PROTEIN')
        
            
    def add_peptide(self, peptide):  
        if 'PROTEIN_ID' in peptide.columns:
            peptide_protein_mapping = peptide[['ID', 'PROTEIN_ID']] \
                .rename(columns={'ID': 'PEPTIDE_ID'})     
            peptide = peptide.drop(columns=['PROTEIN_ID'])
        else:
            peptide_protein_mapping = None
            
        if 'GENE_ID' in peptide.columns:
            peptide_gene_mapping = peptide[['ID', 'GENE_ID']] \
                .rename(columns={'ID': 'PEPTIDE_ID'}) 
            peptide = peptide.drop(columns=['GENE_ID'])
        else:
            peptide_gene_mapping = None
            
        self.add_table(peptide, 'PEPTIDE')
        if peptide_protein_mapping is not None:
            self.add_table(peptide_protein_mapping, 'PEPTIDE_PROTEIN_MAPPING')
        if peptide_gene_mapping is not None:    
            self.add_table(peptide_gene_mapping, 'PEPTIDE_GENE_MAPPING')
                
            
    def add_precursor(self, precursor):                
        if 'PEPTIDE_ID' in precursor.columns:
            precursor_peptide_mapping = precursor[['ID', 'PEPTIDE_ID']] \
                .rename(columns={'ID': 'PRECURSOR_ID'})   
            precursor = precursor.drop(columns=['PEPTIDE_ID'])
        else:
            precursor_peptide_mapping = None            
            
        if 'COMPOUND_ID' in precursor.columns:
            precursor_compound_mapping = precursor[['ID', 'COMPOUND_ID']] \
                .rename(columns={'ID': 'PRECURSOR_ID'})     
            precursor = precursor.drop(columns=['COMPOUND_ID'])
        else:
            precursor_compound_mapping = None
                        
        self.add_table(precursor, 'PRECURSOR')
        if precursor_peptide_mapping is not None:
            self.add_table(
                precursor_peptide_mapping,
                'PRECURSOR_PEPTIDE_MAPPING'
            )
        if precursor_compound_mapping is not None:
            self.add_table(
                precursor_compound_mapping, 
                'PRECURSOR_COMPOUND_MAPPING'
            )
            
        
    def add_transition(self, transition):  
        if 'PRECURSOR_ID' in transition.columns:
            transition_precursor_mapping = transition[['ID', 'PRECURSOR_ID']] \
                .rename(columns={'ID': 'TRANSITION_ID'})    
            transition = transition.drop(columns=['PRECURSOR_ID'])            
        else:
            transition_precursor_mapping = None
            
        self.add_table(transition, 'TRANSITION')
        if transition_precursor_mapping is not None:
            self.add_table(
                transition_precursor_mapping,
                'TRANSITION_PRECURSOR_MAPPING'
            )
       
        
    def prepare_tables(self, data):
        cursor = self.connection.cursor()
        
        protein = data[['ProteinId', 'Decoy']] \
            .drop_duplicates(['ProteinId'])
            
        cursor.execute(sqlconst.SELECT_MAX_ID_PROTEIN)
        protein_max_id = cursor.fetchone()[0]   
        
        if protein_max_id is None:
            protein['Exist'] = False
            protein['ID'] = list(range(0, len(protein)))
            protein_max_id = len(protein)
            
        else:
            protein_id = [
                next((x for i, x in self.get_protein(accession=x)['ID'] \
                      .iteritems()), None)
                for x in protein['ProteinId'].values
            ]
            protein['Exist'] = list(map(lambda x: x is not None, protein_id))
            for i, x in enumerate(protein_id):
                if x is None:
                    protein_max_id += 1
                    protein_id[i] = protein_max_id
                    
            protein['ID'] = protein_id
        
        protein.set_index('ProteinId', drop=False, inplace=True)
        
        
        peptide = data[[
            'ProteinId', 
            'PeptideSequence', 'ModifiedPeptideSequence', 
            'Decoy'
        ]].drop_duplicates(['ModifiedPeptideSequence'])
        
        cursor.execute(sqlconst.SELECT_MAX_ID_PEPTIDE)
        peptide_max_id = cursor.fetchone()[0]
        
        if peptide_max_id is None:
            peptide['Exist'] = False
            peptide['ID'] = list(range(0, len(peptide)))
            peptide_max_id = len(peptide)
        
        else:
            peptide_id = [
                next((x for i, x in self.get_peptide(modified_sequence=x)['ID'] \
                      .iteritems()), None)
                for x in peptide['ModifiedPeptideSequence'].values
            ]
            peptide['Exist'] = list(map(lambda x: x is not None, peptide_id))
            for i, x in enumerate(peptide_id):
                if x is None:
                    peptide_max_id += 1
                    peptide_id[i] = peptide_max_id
                    
            peptide['ID'] = peptide_id
        
        peptide.set_index('ModifiedPeptideSequence', drop=False, inplace=True)
        
                        
        precursor = data[[
            'ModifiedPeptideSequence', 
            'TransitionGroupId', 'PrecursorMz', 'PrecursorCharge',
            'NormalizedRetentionTime', 
            'Decoy'
        ]].drop_duplicates(['TransitionGroupId'])
        
        cursor.execute(sqlconst.SELECT_MAX_ID_PRECURSOR)
        precursor_max_id = cursor.fetchone()[0]
        
        if precursor_max_id is None:
            precursor['Exist'] = False
            precursor['ID'] = list(range(0, len(precursor)))
            precursor_max_id = len(precursor)
        
        else:
            precursor_id = [
                next((x for i, x in self.get_precursor(precursor_name=x)['ID'] \
                      .iteritems()), None)
                for x in precursor['TransitionGroupId'].values
            ]
            precursor['Exist'] = list(map(lambda x: x is not None, precursor_id))
            for i, x in enumerate(precursor_id):
                if x is None:
                    precursor_max_id += 1
                    precursor_id[i] = precursor_max_id
                    
            precursor['ID'] = precursor_id
            
        precursor.set_index('TransitionGroupId', drop=False, inplace=True)
        
                
        transition = data[[
            'TransitionGroupId', 
            'TransitionId', 'ProductMz', 'ProductCharge', 'FragmentType',
            'Annotation', 'FragmentSeriesNumber', 'LibraryIntensity', 
            'DecoyTransition' \
                if 'DecoyTransition' in data.columns \
                else 'Decoy',
            'DetectingTransition', 'IdentifyingTransition',
            'QuantifyingTransition'
        ]].drop_duplicates(['TransitionId'])
        if 'DecoyTransition' in transition.columns:
            transition.rename(
                columns={'DecoyTransition': 'Decoy'}, 
                inplace=True
            )
        
        cursor.execute(sqlconst.SELECT_MAX_ID_TRANSITION)
        transition_max_id = cursor.fetchone()[0]
        
        if transition_max_id is not None:
            for x in data['TransitionId'].values:
                if len(self.get_transition(transition_name=x)) > 0:
                    raise ValueError('transition exists: ' + str(x))
                    
            transition['ID'] = list(range(transition_max_id + 1, \
                transition_max_id + 1 + len(transition)))
            transition_max_id += len(transition)
            
        else:
            transition['ID'] = list(range(0, len(transition)))
            transition_max_id = len(transition)
   
        transition.set_index('TransitionId', drop=False, inplace=True)
        
        gene = self.get_gene(gene_name='NA')
        if len(gene) == 0:
            cursor.execute(sqlconst.SELECT_MAX_ID_PEPTIDE)
            gene_max_id = cursor.fetchone()[0]
            if gene_max_id is None:
                gene_id = 0
            else:
                gene_id = gene_max_id + 1
            gene = pd.DataFrame.from_records(
                [(gene_id, 'NA', 1)],
                columns=['ID', 'GENE_NAME', 'DECOY']
            )       
        else:
            gene_id = gene.loc[0, 'ID']
            gene = None
        
        transition = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['TransitionId'], 
                 x['ProductMz'], x['ProductCharge'], x['FragmentType'],
                 x['Annotation'], x['FragmentSeriesNumber'], 
                 x['DetectingTransition'],
                 x['IdentifyingTransition'],
                 x['QuantifyingTransition'],
                 x['LibraryIntensity'], x['Decoy'],
                 precursor.loc[x['TransitionGroupId'], 'ID'])
                for i, x in transition.iterrows()
            ],
            columns=[
                'ID',
                'TRAML_ID',
                'PRODUCT_MZ',
                'CHARGE',
                'TYPE',
                'ANNOTATION',
                'ORDINAL',
                'DETECTING',
                'IDENTIFYING',
                'QUANTIFYING',
                'LIBRARY_INTENSITY',
                'DECOY',
                'PRECURSOR_ID'
            ]
        )
            
        precursor = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['TransitionGroupId'], 
                 None, # GROUP_LABEL
                 x['PrecursorMz'], x['PrecursorCharge'],
                 None, # LIBRARY_INTENSITY
                 x['NormalizedRetentionTime'], 
                 None, # LIBRARY_DRIFT_TIME
                 x['Decoy'],
                 peptide.loc[x['ModifiedPeptideSequence'], 'ID'])
                for i, x in precursor.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'TRAML_ID',
                'GROUP_LABEL',
                'PRECURSOR_MZ',
                'CHARGE',
                'LIBRARY_INTENSITY',
                'LIBRARY_RT',
                'LIBRARY_DRIFT_TIME',
                'DECOY',
                'PEPTIDE_ID'
            ]
        )
            
        peptide = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['PeptideSequence'], 
                 x['ModifiedPeptideSequence'], x['Decoy'],
                 protein.loc[x['ProteinId'], 'ID'],
                 gene_id)
                for i, x in peptide.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'UNMODIFIED_SEQUENCE',
                'MODIFIED_SEQUENCE',
                'DECOY',
                'PROTEIN_ID',
                'GENE_ID'
            ]   
        )
            
        protein = pd.DataFrame.from_records(
            [
                (int(x['ID']), x['ProteinId'], x['Decoy'])
                for i, x in protein.iterrows()
                if not x['Exist']
            ],
            columns=[
                'ID',
                'PROTEIN_ACCESSION',
                'DECOY'
            ]
        )
        
        return {
            'gene': gene,
            'protein': protein,
            'peptide': peptide,
            'precursor': precursor,
            'transition': transition
        }
        

    def add_data(self, data):
        tables = self.prepare_tables(data)
          
        for k, table in tables.items():
            add_func = getattr(self, 'add_' + k, None)
            if callable(add_func):
                add_func(table)
            else:
                self.add_table(table, k.upper())
        
            
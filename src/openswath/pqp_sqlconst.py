CREATE_TABLE_COMPOUND = '''
CREATE TABLE COMPOUND (
    ID            INT  PRIMARY KEY
                       NOT NULL,
    COMPOUND_NAME TEXT NOT NULL,
    SUM_FORMULA   TEXT NOT NULL,
    SMILES        TEXT NOT NULL,
    ADDUCTS       TEXT NOT NULL,
    DECOY         INT  NOT NULL
);
'''

CREATE_TABLE_GENE = '''
CREATE TABLE GENE (
    ID        INT  PRIMARY KEY
                   NOT NULL,
    GENE_NAME TEXT NOT NULL,
    DECOY     INT  NOT NULL
);
'''

CREATE_TABLE_PEPTIDE = '''
CREATE TABLE PEPTIDE (
    ID                  INT  PRIMARY KEY
                             NOT NULL,
    UNMODIFIED_SEQUENCE TEXT NOT NULL,
    MODIFIED_SEQUENCE   TEXT NOT NULL,
    DECOY               INT  NOT NULL
);
'''

CREATE_TABLE_PEPTIDE_GENE_MAPPING = '''
CREATE TABLE PEPTIDE_GENE_MAPPING (
    PEPTIDE_ID INT NOT NULL,
    GENE_ID    INT NOT NULL
);
'''

CREATE_TABLE_PEPTIDE_PROTEIN_MAPPING = '''
CREATE TABLE PEPTIDE_PROTEIN_MAPPING (
    PEPTIDE_ID INT NOT NULL,
    PROTEIN_ID INT NOT NULL
);
'''

CREATE_TABLE_PRECURSOR = '''
CREATE TABLE PRECURSOR (
    ID                 INT  PRIMARY KEY
                            NOT NULL,
    TRAML_ID           TEXT,
    GROUP_LABEL        TEXT,
    PRECURSOR_MZ       REAL NOT NULL,
    CHARGE             INT,
    LIBRARY_INTENSITY  REAL,
    LIBRARY_RT         REAL,
    LIBRARY_DRIFT_TIME REAL,
    DECOY              INT  NOT NULL
);
'''

CREATE_TABLE_PRECURSOR_COMPOUND_MAPPING = '''
CREATE TABLE PRECURSOR_COMPOUND_MAPPING (
    PRECURSOR_ID INT NOT NULL,
    COMPOUND_ID  INT NOT NULL
);
'''

CREATE_TABLE_PRECURSOR_PEPTIDE_MAPPING = '''
CREATE TABLE PRECURSOR_PEPTIDE_MAPPING (
    PRECURSOR_ID INT NOT NULL,
    PEPTIDE_ID   INT NOT NULL
);
'''

CREATE_TABLE_PROTEIN = '''
CREATE TABLE PROTEIN (
    ID                INT  PRIMARY KEY
                           NOT NULL,
    PROTEIN_ACCESSION TEXT NOT NULL,
    DECOY             INT  NOT NULL
);
'''

CREATE_TABLE_TRANSITION = '''
CREATE TABLE TRANSITION (
    ID                INT      PRIMARY KEY
                               NOT NULL,
    TRAML_ID          TEXT,
    PRODUCT_MZ        REAL     NOT NULL,
    CHARGE            INT,
    TYPE              CHAR (1),
    ANNOTATION        TEXT,
    ORDINAL           INT,
    DETECTING         INT      NOT NULL,
    IDENTIFYING       INT      NOT NULL,
    QUANTIFYING       INT      NOT NULL,
    LIBRARY_INTENSITY REAL,
    DECOY             INT      NOT NULL
);
'''

CREATE_TABLE_TRANSITION_PEPTIDE_MAPPING = '''
CREATE TABLE TRANSITION_PEPTIDE_MAPPING (
    TRANSITION_ID INT NOT NULL,
    PEPTIDE_ID    INT NOT NULL
);
'''

CREATE_TABLE_TRANSITION_PRECURSOR_MAPPING = '''
CREATE TABLE TRANSITION_PRECURSOR_MAPPING (
    TRANSITION_ID INT NOT NULL,
    PRECURSOR_ID  INT NOT NULL
);
'''

CREATE_TABLE_VERSION = '''
CREATE TABLE VERSION (
    ID INT NOT NULL
);
'''


INSERT_VERSION = '''INSERT INTO VERSION (ID) VALUES (?);'''

INSERT_GENE = '''
INSERT INTO GENE (
    ID,
    GENE_NAME,
    DECOY
) VALUES (?,?,?);
'''

INSERT_PEPTIDE = '''
INSERT INTO PEPTIDE (
    ID,
    UNMODIFIED_SEQUENCE,
    MODIFIED_SEQUENCE,
    DECOY
) VALUES (?,?,?,?);
'''

INSERT_PEPTIDE_GENE_MAPPING = '''
INSERT INTO PEPTIDE_GENE_MAPPING (
    PEPTIDE_ID,
    GENE_ID 
) VALUES (?,?);
'''

INSERT_PEPTIDE_PROTEIN_MAPPING = '''
INSERT INTO PEPTIDE_PROTEIN_MAPPING (
    PEPTIDE_ID,
    PROTEIN_ID
) VALUES (?,?);
'''

INSERT_PRECURSOR = '''
INSERT INTO PRECURSOR (
    ID,
    TRAML_ID,
    GROUP_LABEL,
    PRECURSOR_MZ,
    CHARGE,
    LIBRARY_INTENSITY,
    LIBRARY_RT,
    LIBRARY_DRIFT_TIME,
    DECOY
) VALUES (?,?,?,?,?,?,?,?,?);
'''

INSERT_PRECURSOR_PEPTIDE_MAPPING = '''
INSERT INTO PRECURSOR_PEPTIDE_MAPPING (
    PRECURSOR_ID,
    PEPTIDE_ID
) VALUES (?,?);
'''

INSERT_PROTEIN = '''
INSERT INTO PROTEIN (
    ID,
    PROTEIN_ACCESSION,
    DECOY
) VALUES (?,?,?);
'''

INSERT_TRANSITION = '''
INSERT INTO TRANSITION (
    ID,
    TRAML_ID,
    PRODUCT_MZ,
    CHARGE,
    TYPE,
    ANNOTATION,
    ORDINAL,
    DETECTING,
    IDENTIFYING,
    QUANTIFYING,
    LIBRARY_INTENSITY,
    DECOY
) VALUES (?,?,?,?,?,?,?,?,?,?,?,?);
'''

INSERT_TRANSITION_PRECURSOR_MAPPING = '''
INSERT INTO TRANSITION_PRECURSOR_MAPPING (
    TRANSITION_ID,
    PRECURSOR_ID
) VALUES (?,?);
'''


SELECT_GENE = '''
SELECT * FROM GENE
'''

SELECT_PROTEIN = '''
SELECT * FROM PROTEIN
'''

SELECT_PEPTIDE = '''
SELECT * FROM PEPTIDE
'''

JOIN_PEPTIDE_PROTEIN = '''
INNER JOIN (
    SELECT PEPTIDE_ID AS ID, PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING
)
USING (ID)
'''

JOIN_PEPTIDE_GENE = '''
INNER JOIN (
    SELECT PEPTIDE_ID AS ID, GENE_ID FROM PEPTIDE_GENE_MAPPING
)
USING (ID)
'''

SELECT_PRECURSOR = '''
SELECT * FROM PRECURSOR
'''

JOIN_PRECURSOR_PEPTIDE = '''
INNER JOIN (
    SELECT PRECURSOR_ID AS ID, PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING
)
USING (ID)
'''

JOIN_PRECURSOR_COMPOUND = '''
INNER JOIN (
    SELECT PRECURSOR_ID AS ID, COMPOUND_ID FROM PRECURSOR_COMPOUND_MAPPING
)
USING (ID)
'''

SELECT_TRANSITION = '''
SELECT * FROM TRANSITION
'''

JOIN_TRANSITION_PRECURSOR = '''
INNER JOIN (
    SELECT TRANSITION_ID AS ID, PRECURSOR_ID FROM TRANSITION_PRECURSOR_MAPPING
)
USING (ID)
'''


SELECT_MAX_ID_GENE = '''
SELECT MAX(ID) FROM GENE
'''

SELECT_MAX_ID_PROTEIN = '''
SELECT MAX(ID) FROM PROTEIN
'''

SELECT_MAX_ID_PEPTIDE = '''
SELECT MAX(ID) FROM PEPTIDE
'''

SELECT_MAX_ID_PRECURSOR = '''
SELECT MAX(ID) FROM PRECURSOR
'''

SELECT_MAX_ID_TRANSITION = '''
SELECT MAX(ID) FROM TRANSITION
'''


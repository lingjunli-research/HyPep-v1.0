import pandas as pd
import collections
import pickle

from sympy import sequence


with open('database_path.pkl', 'rb') as file_e:      
    # Call load method to deserialze
    db_path = pickle.load(file_e)

with open('peaks_path.pkl', 'rb') as file_l:      
    # Call load method to deserialze
    peaks_path = pickle.load(file_l)

with open('motif_path.pkl', 'rb') as file_c:      
    # Call load method to deserialze
    motif_path = pickle.load(file_c)

def input_list():
    """
    ### Retrieve data from user's PEAKS_export.csv input file

    Return:
        Pandas Dataframe with headers: Peptide, Scan, ALC(%), PTM, and Local Confidence(%)
    """
    return pd.read_csv(peaks_path, usecols=['Peptide','Scan','ALC (%)','PTM','local confidence (%)'])


def input_database():
    return pd.read_csv(db_path, usecols=['Accession Number', 'Sequence'])


def database_sequence_data():
    """
    ### Retrieve data from user's database input file
    ### includes L's and I's modification

    Return:
        database_sequence (dict): a dictionary consists of keys as the modified database sequence and values as the original database sequence
    """
    database_df = input_database()
    db_sequence = list(database_df['Sequence'])
    database_sequence = collections.defaultdict(str)
    #database_with_accession = collections.defaultdict(str)
    for db_seq in db_sequence:
        while '(' in db_seq: 
            starting = db_seq.index('(')
            ending = db_seq.index(')')
            ptm = db_seq[starting:ending+1]
            db_seq = db_seq.replace(ptm,'') 
        no_ptm_db = ''
        for amino_acid in db_seq:
            if amino_acid.isupper():
                no_ptm_db+=amino_acid
        database_sequence[no_ptm_db] = no_ptm_db.replace('I','L')
        #accession_number = str(database_df.loc[database_df['Sequence'] == db_seq].iloc[0,0] + ',')
    return database_sequence

def database_accession_data():
    """

    Return:
        database_sequence (dict): a dictionary consists of keys as the modified database sequence and values as the original database sequence
    """
    database_df = input_database()
    db_sequence = list(database_df['Sequence'])
    accession_num = list(database_df['Accession Number'])
    accession_db_dict = collections.defaultdict(str)
    for seq in range(len(db_sequence)):
        db_seq = db_sequence[seq]
        while '(' in db_seq: 
            db_seq = db_seq.replace(db_seq[db_seq.index('('):db_seq.index(')')+1],'') 
        no_ptm_db = ''
        for amino_acid in db_seq:
            if amino_acid.isupper():
                no_ptm_db+=amino_acid
        accession_db_dict[no_ptm_db] += str(str(accession_num[seq])+',')
    return accession_db_dict

#for key, val in database_accession_data().items():
#    print(key, val)

def database_list():
    """
    ### A list of original database sequence
    
    Return:
        List of database sequences without I/L modifications
    """
    return list( database_sequence_data().keys())


def precursor_sequence_data():
    """
    ### Retrieve and extract de novo sequenced peptides data, scan number, and ALC % from user's PEAKS_export.csv input file
    - PTMs are removed from de novo sequenced peptides
    - includes L's and I's modification

    Return:
        modified_precursor_sequence_dict (dict): a dictionary consists of keys as the modified precursor sequence and values as the original precursor sequence
        scan_number_dict (dict): a dictionary consists of keys as precursor sequence and values as their corresponding scan number and ALC percentages
    """
    de_novo_sequenced_peptide_data = input_list()
    peptide_list = list(de_novo_sequenced_peptide_data['Peptide'])
    modified_precursor_sequence_dict = collections.defaultdict(str)
    scan_number_dict = collections.defaultdict(str)
    for i in range(len(peptide_list)):
        sequence = peptide_list[i]
        while '(' in sequence: 
            starting = sequence.index('(')
            ending = sequence.index(')')
            ptm = sequence[starting:ending+1]
            sequence = sequence.replace(ptm,'')
        scan_info = str(de_novo_sequenced_peptide_data.iloc[i,1]) + '[' +  str(de_novo_sequenced_peptide_data.iloc[i,2]) + '],'
        scan_number_dict[sequence] += scan_info
        modified_precursor_sequence_dict[sequence] = sequence.replace('I','L')
    return modified_precursor_sequence_dict, scan_number_dict

def motifs():
    """
    ### Retrieve data from user's motif.csv input file

    Return:
        Pandas Dataframe with headers: motif sequences and the starting positions
    """
    return pd.read_csv(motif_path, usecols=['Motif Sequence', 'Start Position'])
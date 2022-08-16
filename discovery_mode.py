import identification_mode
from input_data import precursor_sequence_data, input_list, motifs

def precursor_sorted(SHS_data=list):
    """
    ### Sorting Precursor Sequence from A-Z and Score from largest to smallest and remove the duplicated precursor sequence

    Args:
        SHS_data (list): SHS raw output
    Return:
        precursor_sequence (list): all non-duplicated precursor sequence in a list (nested list)
    """  
    
    peptide_list = identification_mode.ID_modifications.retrieving_SHS_data(SHS_data)
    pr_sort = sorted(peptide_list, key = lambda x: (x[0],-float(x[3])))
    p_sequence = []
    pr_check = []
    for pr_sequence in pr_sort:
        if float(pr_sequence[3]) >= 4.0:
            continue
        if pr_sequence[0] not in pr_check:
            p_sequence.append(pr_sequence)
            pr_check.append(pr_sequence[0])
    return p_sequence    

def ALC_filter(SHS_data=list, ALC=float):
    """
    ### Filtering the SHS data by user-defined ALC percentage

    Args:
        SHS_data (list): SHS output at user-defined FDR
        ALC (float): user-defined ALC
    Return:
        sorted_list (list): SHS data list with ALC% filtered
    """  
    alc_filtered_list = []
    for data in SHS_data:
        scan_alc_filtered_string= ''
        scan_num_alc = data[1].split(',')
        for scan in scan_num_alc:
            if len(scan)==0:
                continue
            alc_percentage = scan[int(scan.index('[')+1):scan.index(']')]
            scan_number = scan[0:scan.index('[')]
            if float(alc_percentage) >= ALC:
                scan_alc_filtered_string+=str(scan_number+'['+alc_percentage+'],')
        
        
        if len(scan_alc_filtered_string)!=0:
            data[1] = scan_alc_filtered_string[0:-1]
            alc_filtered_list.append(data)
    sorted_list = sorted(alc_filtered_list, key = lambda x: (x[2],-float(x[3])))
    return sorted_list

def local_confidence_check(user_local_confidence, local_confidence_list):
    """
    ### Identify if all local confidence from the precursor sequences is greater than user-defined local confidence

    Args:
        user_local_confidence (list): user-defined local confidence percentage
        local_confidence_list (float): all local confidence percentage with their corresponding precursor sequences
    Return:
        True if all local confidence percentages are greater than user-defined local confidence percentage
        False if any local confidence percentage is less than user-defined local confidence percentage 
    """ 
    for LC in local_confidence_list:
        if user_local_confidence >= int(LC):
            return False
    return True

def local_confidence_filter(SHS_FDR_filtered_data=list, user_local_confidence=float):
    """
    ### Filtering the SHS data by user-defined local confidence percentage

    Args:
        SHS_FDR_filtered_data (list): SHS output at user-defined FDR
        user_local_confidence (float): user-defined local confidence percentage
    Return:
        LC_filtered_list (list): SHS data list with local confidence % filtered
    """
    precursor_dataframe = input_list()
    LC_filtered_list = []
    for data in SHS_FDR_filtered_data:
        scans_info = ''
        scan_list = data[1].split(',')
        scan_pass_list = []
        for scans_alc in scan_list:
            scans = scans_alc[0:scans_alc.index('[')]
            local_confidence_list = precursor_dataframe.loc[precursor_dataframe['Scan'] == scans].iloc[0,3].split()
            if local_confidence_check(user_local_confidence,local_confidence_list):
                scan_pass_list.append(scans_alc)
        for i in scan_pass_list:
            scans_info+=str(i+',')
        if len(scans_info) != 0:
            data[1]= scans_info[0:-1]
            LC_filtered_list.append(data)
    return LC_filtered_list

def PTM(precursor_sequence, ptm_description_list):
    """
    ### Replace the monoisotopic mass in the de novo sequenced peptide to actual PTMs

    Args:
        precursor_sequence (str): de novo sequenced peptide
        ptm_description_list (list): PTMs description
    Return:
        precursor sequence (str): de novo sequenced peptide with monoisotopic mass replaced with modifications
    """
    if ptm_description_list == ['nan']:
        return precursor_sequence
    PTM_dict = {'Amidation': '(-.98)',
                'Oxidation (M)': '(+15.99)',
                'Pyro-glu from Q':'Q(-17.03)',
                'Pyro-glu from E':'E(-18.01)'}
    PTM_replacement = {'Amidation': 'amide',
                        'Oxidation (M)': '(O)',
                        'Pyro-glu from Q':'pQ',
                        'Pyro-glu from E':'pE'}
    
    while len(ptm_description_list) > 0:
        ptm_des = ptm_description_list[0]
        precursor_sequence = precursor_sequence.replace(PTM_dict[ptm_des], PTM_replacement[ptm_des])
        ptm_description_list.remove(ptm_description_list[0])
    return precursor_sequence
    
def motif_indicator(motif=str, position=str, precursor_sequece=str):
    """
    ### Searches the motif sequence in the de novo sequenced peptide
    
    Args:
        motif (str): motif sequence
        position (str): starting search position ('C-term', 'N-term', or 'Any')
        precursor_sequence (str): de novo sequenced peptide
    Return:
        True if motif sequence found in the de novo sequenced peptide
        False if motif sequence is not found in the de novo sequenced peptide
    """
    if position == 'C-term' and precursor_sequece.endswith(motif):
        return True
    if position == 'N-term' and precursor_sequece.startswith(motif):
        return True
    if position == 'Any' and motif in precursor_sequece:
        return True
    else:
        return False

def motif_filter(SHS_FDR_filtered_data):
    """
    ### Reports if motif sequences are found in the precursor sequences 
    
    Args:
        SHS_FDR_filtered_data (list): SHS data at user-defined FDR

    Return:
        motif_filtered_list (list): a list contanins matches with motif sequences found in the precursor sequence
    """ 
    precursor_dataframe = input_list()
    motif = motifs()
    motif_filtered_list = [] 
    motif_sequences = list(motif['Motif Sequence'])
    for data in SHS_FDR_filtered_data:
        scans_info = ''
        scan_list = data[1].split(',')
        for scanning in scan_list:
            scans = scanning[0:scanning.index('[')]
            precursor_sequence, ptm_description = [precursor_dataframe.loc[precursor_dataframe['Scan'] == scans].iloc[0,0], precursor_dataframe.loc[precursor_dataframe['Scan'] == scans].iloc[0,3]] 
            ptm_modified_precursor_sequence = PTM(precursor_sequence, str(ptm_description).split('; '))
            for m in motif_sequences:
                if motif_indicator(m, motif.loc[motif['Motif Sequence'] == m].iloc[0,1], ptm_modified_precursor_sequence):
                    scans_info += str(scanning+',')
        if len(scans_info) != 0:
            data[1] = scans_info[0:-1]
            motif_filtered_list.append(data)
    return motif_filtered_list

            
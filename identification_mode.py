import numpy as np
#from regex import I
from input_data import database_accession_data
import pandas as pd
import pickle

class ID_modifications:
    '''removing duplicated matches from precursor and database sequence'''
    def retrieving_SHS_data(SHS_data):
        """
        ### Receiving data from raw SHS output data and removing the headers and empty spaces

        Args:
            SHS_data (list): SHS raw output
        Return:
            space_removed_peptide_list (list): all SHS data in list with spaced and headers removed
        """

        #peptide_list = [line.split() for line in SHS_data]
        #space_removed_peptide_list = [ele for ele in peptide_list if ele != []]
        #if ['Precursor_Sequence','Scan_Numbers[ALC%]','Database_Sequence','Score'] in space_removed_peptide_list:
        #    space_removed_peptide_list.remove(['Precursor_Sequence','Scan_Numbers[ALC%]','Database_Sequence','Score'])
        #return space_removed_peptide_list
        return SHS_data

    def database_sorted(SHS_data):
        """
        ### Sorting Database Sequence from A-Z and Score from largest to smallest and remove the duplicated database sequence
        
        Args:
            SHS_data (list): SHS raw output
        Return:
            db_sorted_sequence (list): all non-duplicated database sequence in a list (nested list)
        """
        peptide_list = ID_modifications.retrieving_SHS_data(SHS_data) 
        db_sort = sorted(peptide_list, key = lambda x: (x[2],-float(x[3])))
        db_sorted_sequence = []
        db_check = []
        for db_sequence in db_sort:
            if db_sequence[2] not in db_check:
                db_sorted_sequence.append(db_sequence)
                db_check.append(db_sequence[2])
        return db_sorted_sequence
        
    def precursor_sorted(SHS_data):
        """
        ### Sorting Precursor Sequence from A-Z and Score from largest to smallest and remove the duplicated precursor sequence

        Args:
            SHS_data (list): SHS raw output
        Return:
            precursor_sequence (list): all non-duplicated precursor sequence in a list (nested list)
        """ 
        peptide_list = ID_modifications.database_sorted(SHS_data)
        pr_sort = sorted(peptide_list, key = lambda x: (x[0],-float(x[3])))
        precursor_sequence = []
        pr_check = []
        for pr_sequence in pr_sort:
            if pr_sequence[0] not in pr_check:
                precursor_sequence.append(pr_sequence)
                pr_check.append(pr_sequence[0])
        return precursor_sequence
    
    def accesssion_information(SHS_data):
        accession_number = database_accession_data()
        #peptide_list = ID_modifications.retrieving_SHS_data(SHS_data)
        list_with_db_accession = []
        for matches in SHS_data:
            accession_num = accession_number[matches[2]]
            new_accession_information = [matches[0],matches[1], matches[2],accession_num[:-1], matches[3]]
            list_with_db_accession.append(new_accession_information)

        with open('SHS_import_to_AMM.pkl', 'wb') as file_a:
            pickle.dump(list_with_db_accession, file_a)
        return list_with_db_accession   
    def removing_matches(SHS_data, Discovery_data):
        precursor_with_only_one_scan = []
        filtered_discovery_data = []
        for data in SHS_data:
            scan_info = data[1].split(',')
            if len(scan_info) == 1:
                precursor_with_only_one_scan.append(data[0])
        for dis_data in Discovery_data:
            if dis_data[0] in precursor_with_only_one_scan:
                continue
            else:
                filtered_discovery_data.append(dis_data)
        return filtered_discovery_data
                     
class FDR_Calc:
    """FDR Calculation Methods"""
    def FDR_calculation(num_of_target,num_of_decoy):
        """
        %FDR = (decoy_hits/target_hits) * 100
        
        Args:
            num_of_target (int): number of target hits
            num_of_decoy (int): number of decoy hits
        Return:
            Percent FDR (float)
        """

        return round(((num_of_decoy/(num_of_target+num_of_decoy))*100),2)
    
    def filter_by_score(lst, score):
        """
        Filtering the SHS data by algorithm score
        
        Args:
            lst (list): SHS data
            score (float): SHS algorithm score
        Return:
            filtered_list (list): SHS data list filtered by score
        """
        
        filtered_list_update = []

        #list_df = pd.DataFrame([sub.split(" ") for sub in lst])

        #list_df[3] = list_df[3].apply(pd.to_numeric,errors='coerce')
        #list_df_filtered = list_df[list_df[3] >= score]
        
        #precursor_ID = list_df_filtered[0].values.tolist()
        #scan_ID = list_df_filtered[1].values.tolist()
        #db_ID = list_df_filtered[2].values.tolist()
        #score_ID = list_df_filtered[3].values.tolist()
        
        #top_range = len(list_df_filtered)
        
        #for a in range(0,top_range):
            
            #nest = []
            
            #precursor_ID_select = precursor_ID[a]
            #scan_ID_select = scan_ID[a]
            #db_ID_select = db_ID[a]
            #score_ID_select = score_ID[a]
            
            #nest.append(precursor_ID_select)
            #nest.append(scan_ID_select)
            #nest.append(db_ID_select)
            #nest.append(score_ID_select)
            
            #filtered_list.append(nest)
        for data in lst:
            if float(data[3]) >= score:
                filtered_list_update.append(data) 
        return filtered_list_update
        
    def find_closest_FDR(fdr_list, user_fdr):
        """
        Finding the closest actual FDR from all possible FDRs.
        
        Args:
            fdr_list (list): all possible FDR 
        """
        np_fdr_list = np.array(fdr_list)
        difference = np.absolute(np_fdr_list-user_fdr)
        fdr_index=difference.argmin()
        return fdr_index, np_fdr_list[fdr_index]

    def FDR_filter(user_fdr, target_list, decoy_list):
        """
        Removing theoretical false positives by user FDR threshold

        Args:
            user_fdr (float): user-defined FDR
            target_list (list): raw SHS target data
            decoy_list (list): raw SHS decoy data
        Return:
            closest_fdr (float): the closest actual FDR with the user-defined FDR
            fdr_score (float): algorithm score corresponds with the actual FDR
            FDR_filtered_target (list): target list filtered by FDR and its corresponding algorithm score
        """
        decoy_score_list = []
        for decoys in decoy_list:
            decoy_score_list.append(float(decoys[3]))
        decoy_score_list = sorted(list(set(decoy_score_list)),reverse=True)
        FDR_list = []
        for decoy_score in decoy_score_list:
            filter_target_score = FDR_Calc.filter_by_score(target_list,decoy_score)
            filter_decoy_score = FDR_Calc.filter_by_score(decoy_list,decoy_score)
            fdr_percentage = FDR_Calc.FDR_calculation(len(filter_target_score),len(filter_decoy_score))
            if fdr_percentage >5.0: ###
                break
            FDR_list.append(fdr_percentage) 
        score_index, closest_fdr = FDR_Calc.find_closest_FDR(FDR_list,user_fdr)
        fdr_score = decoy_score_list[score_index]
        FDR_filtered_target = FDR_Calc.filter_by_score(target_list,fdr_score)
        return closest_fdr, fdr_score, FDR_filtered_target

    def FDR_list(target_list, decoy_list):
        """
        Reports all possible FDRs from the decoy data
        
        Args:
            target_list (list): SHS target data
            decoy_list (list): SHS decoy data
        Return:
            FDR_list (list): list contains information of each possible FDR ==> [number of target-hits, number of decoy-hits, FDR percentage, algorithm score corresponded to FDR percentage]
        """
        decoy_score_list = []
        for decoys in decoy_list:
            decoy_score_list.append(float(decoys[3]))
        decoy_score_list = sorted(list(set(decoy_score_list)),reverse=True)
        FDR_list = []
        for decoy_score in decoy_score_list:
            filter_target_score = FDR_Calc.filter_by_score(target_list,decoy_score)
            filter_decoy_score = FDR_Calc.filter_by_score(decoy_list,decoy_score)
            fdr_percentage = FDR_Calc.FDR_calculation(len(filter_target_score),len(filter_decoy_score))
            FDR_list.append([len(filter_target_score),len(filter_decoy_score),fdr_percentage, decoy_score]) 
        return FDR_list
  


    


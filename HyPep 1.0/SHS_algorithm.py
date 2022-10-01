import itertools
from input_data import precursor_sequence_data


class algorithm:
    def sequence_length(precursor,database):
        """
        ### Reports the average length of precursor sequence and the database sequence
        """
        return (len(precursor) + len(database))/2
    
    def fixed_algorithm(precursor_database_compare):
        """
        ### SHS Forward and Backward-Fixed Algorithm
        
        Args:
            precursor_database_compare (list): zipped objects of each amino residues from de novo sequenced peptide and database sequence
        Return:
            alignment score (float): algorithm score reflecting local alignment matches between de novo sequenced peptide and database sequence
        """
        alignment_score = 0
        for precursor, database in precursor_database_compare:
            if precursor == database:
                alignment_score+=1                      
        return alignment_score
    
    def varying_algorithm(precursor, database, varying_site):
        """
        ### SHS Forward and Backward-Varying Algorithm
        
        Args:
            precursor (str): de novo sequenced peptide
            database (str): database sequence
            varying_site (int): user_defined sliding window size
        Return:
            alignment score (float): algorithm score reflecting local alignment matches between de novo sequenced peptide and database sequence
        """
        matched_scanning = 0
        alignment_score = 0
        for i in range(len(precursor)):
            current_precursor_amino_acid = precursor[i]
            scanning_site = database[matched_scanning:(matched_scanning+int(varying_site)+1)]
            if current_precursor_amino_acid in scanning_site:
                forward_occurrance_index = scanning_site.index(current_precursor_amino_acid)
                matched_scanning += forward_occurrance_index + 1
                alignment_score += 1
            else:
                if i == matched_scanning:
                    matched_scanning+=1
                else:
                    matched_scanning+=0
        return alignment_score

    def main_algorithm(sws, database_dict):
        """
        ### SHS Main Algorithm
        
        Args:
            sws (int): user-defined sliding window size
            database_dict (dict): a dictionary consists of keys as the modified database sequence and values as the original database sequence
        Return:
            output data (list): SHS output list contains: [precursor sequence, scan numbers, database sequence, algorithm score]
        """ 
        output_data = []
        
        modified_precursor_dict, scan_number_dict = precursor_sequence_data()
        precursor_data = sorted(modified_precursor_dict.keys())
        
        modified_database_dict = database_dict
        database_data = list(modified_database_dict.keys())
         
        precursor_database_combination = list(itertools.product(precursor_data, database_data))
        
        for original_precursor,original_database in precursor_database_combination:
            
            precursor = modified_precursor_dict[original_precursor]
            database = modified_database_dict[original_database]
          
            algorithm_score = 0

            reverse_precursor = precursor[::-1]
            reverse_database = database[::-1] 
            forward_fixed_precursor_database_compare = list(zip(precursor,database))
            backward_fixed_precursor_database_compare = list(zip(reverse_precursor,reverse_database))

            f_fixed_match_count = algorithm.fixed_algorithm(forward_fixed_precursor_database_compare)
            b_fixed_match_count = algorithm.fixed_algorithm(backward_fixed_precursor_database_compare)
            
            f_varying_match_count = algorithm.varying_algorithm(precursor,database,sws)
            b_varying_match_count = algorithm.varying_algorithm(reverse_precursor,reverse_database,sws) 

            algorithm_score = (f_fixed_match_count+b_fixed_match_count+f_varying_match_count+b_varying_match_count)/algorithm.sequence_length(precursor,database)

            if algorithm_score >= 1:
                scan_num_alcs = scan_number_dict[original_precursor][0:-1]
                output_data.append(str(original_precursor + ' '+scan_num_alcs + ' ' + original_database + ' '+ str(round(algorithm_score,3))))
            
        return output_data
   



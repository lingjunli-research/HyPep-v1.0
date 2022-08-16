import random
import collections 
import pandas as pd
from input_data import database_list

database = database_list() 
class properties:
    def sequence_check(sequence):
        """
        ### Remove all lower case amino acid or PTMs.
        Args:
            sequence (str): original peptide sequence
        Return:
            new_sequence (str): sequnce with lower cases removed 
        """
        if sequence.isupper():
            return sequence
        else:
            new_sequence = ''
            for amino_acid in sequence:
                if amino_acid.isupper():
                    new_sequence+=amino_acid
            return new_sequence
        
    def peptide_length_frequency():
        """
        ### Get all the database sequence length
        
        Return:
            peptide_length_freq_dict (dict): a dictionary with keys as peptide lengths and values as the frequency of the peptide length
        """
        peptide_length_freq_dict = collections.defaultdict(int)
        for database_sequence in database:
            peptide_length_freq_dict[len(properties.sequence_check(database_sequence))]+=1
        return peptide_length_freq_dict

    def amino_acid_list():
        """
        ### Get all the amino acids from the database sequence

        Return: 
            amino_acid_whole_list (list): a list contains all the amino acid from the database sequences
        """
        amino_acid_whole_list = []
        for database_sequence in database:
            database_sequence = properties.sequence_check(database_sequence)
            for amino_acid in database_sequence:
                amino_acid_whole_list.append(amino_acid)
        return amino_acid_whole_list
    
    def all_amino_acids():
        """
        Return:
            list: a list of amino acid
        """
        return list(set(properties.amino_acid_list()))
    
    
class Identity_Threshold:
    def ID_threshold_check(decoy_seqeunce):
        """
        ### Local alignment with the decoy sequence and all possible target sequences in the database.
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            True: if the decoy sequence has less than 50% similarity with all target sequences in the database
            False: if the decoy sequence has more than 50% similarity with any of target sequence in the database
        """
        filtered_target = []
        for t in database:
            if len(decoy_seqeunce) == len(t):
                filtered_target.append(t)
        for target_sequence in filtered_target:
            alignment_score = 0
            for amino_acid in list(zip(target_sequence,decoy_seqeunce)):
                if amino_acid[0] == amino_acid[1]:
                    alignment_score+=1
            if alignment_score/len(decoy_seqeunce) < 0.5:
                continue
            else:
                return False
        return True
    
    def max_ID_threshold(decoy_sequence):
        """
        ### Find the maximum identity threshold percentage from the decoy sequence taht matches to target sequences
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            list:   max_ID_T (float): maximum identity threshold percentage
                    target_index (int): index of corresponding target sequence that yields the maximum identity threshold percentage
        """
        max_ID_T = 0
        target_list_all = database
        target_list = []
        for t in target_list_all:
            if len(t) == len(decoy_sequence):
                target_list.append(t)
        for target_index in range(len(target_list)):
            target_sequence = target_list[target_index]
            alignment_score = 0
            for amino_acid in list(zip(target_sequence,decoy_sequence)):
                if amino_acid[0] == amino_acid[1]:
                    alignment_score+=1
            alignment_ID_threshold = alignment_score/len(decoy_sequence)
            if alignment_ID_threshold > max_ID_T:
                max_ID_T = alignment_ID_threshold
        return [max_ID_T, target_index]
    
    def alignment_sequence_index(decoy_sequence):
        """
        ### Gives the matched residues index from target and decoy sequences
        
        Args:
            decoy_sequence (str): decoy sequence generated from the target-decoy method
        Return:
            match_index_list (list): list contains the index of matched residues from both target and decoy sequences
        """
        max_num, index = Identity_Threshold.max_ID_threshold(decoy_sequence)
        target_sequence = database[index]
        x = 0
        match_index_list = []
        for amino_acid in list(zip(target_sequence,decoy_sequence)):
            if amino_acid[0] == amino_acid[1]:
                match_index_list.append(x)
            x+=1
        return match_index_list 
            
class decoy_database:
    '''Generation of decoy list based off user-defined target-decoy methods'''  
    def reversed_decoy_database():
        """
        ### Target sequences are reversed and generated as the decoy sequences
        
        Return:
            reversed_list (list): Reverse decoy database 
        """
        reversed_list = []
        for database_sequence in database:
            reversed_list.append(properties.sequence_check(database_sequence[::-1]))   
        return reversed_list

    def shuffled_decoy_database():
        """
        ### Residues in each target sequences are randomly shuffled and generated as the decoy sequences
        
        Return:
            shuffled_decoy_list (list): Shuffle decoy database with Identify Threshold applied
        """
        shuffled_decoy_list = []
        
        def shuffled_algorithm(database_sequence):
            """
            ### Shuffles each amino acid residue in each target sequences
            
            Args:
                database_seqeuence (str): target sequence from the database
            
            Return:
                shuffled_database_sequence (str): decoy seqeunce from shuffle decoy method
            """
            shuffled_database_sequence = ''
            for i in random.sample(range(0,len(database_sequence)),len(database_sequence)):
                shuffled_database_sequence += database_sequence[i]
            return shuffled_database_sequence

        shuffled_table = []
        # original shuffled decoy list
        for database_sequence in database:
            shuffled_decoy_sequence = shuffled_algorithm(database_sequence)
            shuffled_decoy_list.append(shuffled_decoy_sequence)
            shuffled_table.append([shuffled_decoy_sequence, round(Identity_Threshold.max_ID_threshold(shuffled_decoy_sequence)[0],3)])
            
        iterations = 1
        over_threshold_list = []
        for decoy_index in range(len(shuffled_decoy_list)):
            decoy_sequence = shuffled_decoy_list[decoy_index]
            if len(decoy_sequence) == 2:    #bypassing sequence with length of 2
                continue
            if not Identity_Threshold.ID_threshold_check(decoy_sequence):
                over_threshold_list.append([decoy_index, decoy_sequence])
        
        while iterations < 31:
            if len(over_threshold_list) == 0:
                break
            reshuffled_list = []
            for index,sequence in over_threshold_list:   
                reshuffled = shuffled_algorithm(sequence)
                reshuffled_list.append([index, sequence, reshuffled])
            for index,sequence,reshuffled in reshuffled_list:
                if Identity_Threshold.ID_threshold_check(reshuffled):
                    shuffled_decoy_list[index] = reshuffled
                    over_threshold_list.remove([index,sequence])
            if len(over_threshold_list) == 0:
                iterations = 30
            iterations+=1
            
        if len(over_threshold_list) > 0:
            all_amino_acid_left = []
            for i, seq in over_threshold_list:
                for s in seq:
                    all_amino_acid_left.append(s)
        
            if len(over_threshold_list) > 0:
                for index, sequence in over_threshold_list:
                    while not Identity_Threshold.ID_threshold_check(sequence):
                        index_list = Identity_Threshold.alignment_sequence_index(sequence)
                        mutation = sequence.replace(sequence[index_list[0]],random.choice(all_amino_acid_left))        ### double check this
                        if Identity_Threshold.ID_threshold_check(mutation):
                            sequence = mutation
                            shuffled_decoy_list[index] = mutation
                        else:
                            continue
        return shuffled_decoy_list

    def randomized_decoy_database():
        """
        ### Amino acids residues are randomly generated to form the decoy sequence
        
        Return:
            randomized_list (list): Random decoy database with Identify Threshold applied
        """
        randomized_list = []

        def randomized_algorithm(all, pep_length): 
            """
            ### Amino acids are randomly cut off from the all amino acids existed in the database. Random decoy sequences are generated based on the frequency of amino acids and peptide length in the database
            
            Args:
                all (list): list that contains the remaining amino acids residues in the database
                pep_length (int): peptide length for generating the decoy sequence
            
            Return:
                sequence (str): decoy seqeunce from random decoy method
            """
            sequence= ''.join(random.sample(all, int(pep_length)))
            return sequence
        

        peptide_length_frequencey = properties.peptide_length_frequency()
        all_database_sequence = properties.amino_acid_list()
        for peptide_length in sorted(peptide_length_frequencey.keys()):
            peptide_length_freq = peptide_length_frequencey[peptide_length] 
            for i in range(peptide_length_freq):
                randomized_decoy_sequence = randomized_algorithm(all_database_sequence,peptide_length)
                for r in randomized_decoy_sequence:
                    all_database_sequence.remove(r)
                randomized_list.append(randomized_decoy_sequence)
        
        iterations = 1
        
        over_threshold_list = []
        overthreshold_all_sequence = []
        for decoy_index in range(len(randomized_list)):
            decoy_sequence = randomized_list[decoy_index]
            if len(decoy_sequence) == 2:    #bypassing sequence with length of 2
                continue
            if Identity_Threshold.ID_threshold_check(decoy_sequence):
                continue
            else:
                over_threshold_list.append([decoy_index, decoy_sequence])
                for d in decoy_sequence:
                    overthreshold_all_sequence.append(d)
        
        if len(over_threshold_list) > 0:#check if there is any decoy sequences that are over 50% identity threshold
            while iterations < 31:
                total=0
                re_random_list = []
                for index,sequence in over_threshold_list:
                    re_random = randomized_algorithm(overthreshold_all_sequence,len(sequence))
                    total+=1
                    for r in re_random:
                        overthreshold_all_sequence.remove(r)
                    re_random_list.append([index,sequence,re_random])
                    
                for index,sequence, re_random in re_random_list:
                    if Identity_Threshold.ID_threshold_check(re_random):
                        randomized_list[index] = re_random
                        over_threshold_list.remove([index,sequence])
                    else:
                        for ran in re_random:
                            overthreshold_all_sequence.append(ran)
    
                if len(over_threshold_list) == 0:
                    iterations = 30
                    break
                iterations+=1
                
        if len(over_threshold_list) > 0:
            all_amino_acid_left = []
            for i, seq in over_threshold_list:
                for s in seq:
                    all_amino_acid_left.append(s)
                    
            for index, sequence in over_threshold_list:
                while not Identity_Threshold.ID_threshold_check(sequence):
                    index_list = Identity_Threshold.alignment_sequence_index(sequence)
                    mutation = sequence.replace(sequence[index_list[0]],random.choice(all_amino_acid_left))
                    if Identity_Threshold.ID_threshold_check(mutation):
                        sequence = mutation
                        randomized_list[index] = mutation
                    else:
                        continue 
        return randomized_list
    
    def hybrid_decoy_database():
        """
        ### Decoy sequences are genereated based on the amino acid residue pattern from the target sequence 
        
        Return:
            randomized_list (list): Random decoy database with Identify Threshold applied 
        """
        decoy_list = []
        all_amino_acid = properties.amino_acid_list() 
        all_amino_acid = properties.amino_acid_list()
        all_database_sequence = sorted(database,key=len)
        
        def hybrid_algorithm(target_sequence=str, all_amino_acid=list):
            """
            ### Amino acids are generated randomly but follow the residue pattern from the target sequence
            
            Args:
                target_sequence (str): sequences from the database
                all_amino_acid (list): all the remaining amino acids from the database
            
            Return:
                sequence (str): decoy seqeunce from hybrid decoy method
            """
            index_list = []
            sequence_index = []
            sequence_index_dict = collections.defaultdict(str)
            for i in target_sequence:
                if i not in index_list:
                    index_list.append(i)
                index = index_list.index(i)
                sequence_index_dict[index] = random.choice(all_amino_acid)
                sequence_index.append(index)
            sequence = ''
            for each_index in sequence_index:
                residue = sequence_index_dict[each_index]
                if residue in all_amino_acid:
                    sequence+=residue
                    all_amino_acid.remove(residue)
                else:
                    sequence_index_dict[each_index] = random.choice(all_amino_acid)
                    sequence+=sequence_index_dict[each_index]
                    all_amino_acid.remove(sequence_index_dict[each_index])
            return sequence, all_amino_acid
        
        decoy_table = []
        for t_sequence in all_database_sequence:
            decoy_seq, all_amino_acid = hybrid_algorithm(t_sequence, all_amino_acid)
            decoy_list.append(decoy_seq)
            decoy_table.append([decoy_seq,round(Identity_Threshold.max_ID_threshold(decoy_seq)[0],3)])

        over_threshold_list = []
        overthreshold_residue_list = []

        for decoy_index in range(len(decoy_list)):
            decoy_seq = decoy_list[decoy_index]
            if len(decoy_seq) == 2:    #bypassing sequence with length of 2
                continue
            if not Identity_Threshold.ID_threshold_check(decoy_seq):
                over_threshold_list.append([decoy_index, decoy_seq])
                for ds in decoy_seq:
                    overthreshold_residue_list.append(ds)
                    
        iterations = 1
        if len(over_threshold_list) > 0:#check if there is any decoy sequences that are over 50% identity threshold
            while iterations < 31:
                re_decoy_list = []
                re_decoy_table = []
                passed = 0 #####
                total = 0
                for index,sequence in over_threshold_list:
                    re_run, overthreshold_residue_list = hybrid_algorithm(all_database_sequence[index],overthreshold_residue_list)
                    re_decoy_list.append([index,sequence, re_run])
                    total += 1
                    re_decoy_table.append([total,re_run,round(Identity_Threshold.max_ID_threshold(re_run)[0],3), Identity_Threshold.ID_threshold_check(re_run)])
                    
                for index,sequence,re_run in re_decoy_list:
                    if Identity_Threshold.ID_threshold_check(re_run):
                        decoy_list[index] = re_run
                        over_threshold_list.remove([index,sequence])
                    else:
                        passed+=1
                        for r in re_run:
                            overthreshold_residue_list.append(r)
                            
                if len(over_threshold_list) == 0:
                    iterations = 30
                    break
                iterations+=1
                
        if len(over_threshold_list) > 0:
            all_amino_acid_left = []
            for i, seq in over_threshold_list:
                for s in seq:
                    all_amino_acid_left.append(s)
                    
            for index, sequence in over_threshold_list:
                while not Identity_Threshold.ID_threshold_check(sequence):
                    index_list = Identity_Threshold.alignment_sequence_index(sequence)
                    mutation = sequence.replace(sequence[index_list[0]],random.choice(all_amino_acid_left))
                    if Identity_Threshold.ID_threshold_check(mutation):
                        sequence = mutation
                        decoy_list[index] = mutation
                    else:
                        continue 
        return decoy_list
    
    def default():
        return 'Incorrect Decoy Database Type'
    
    decoy_database_switcher = {
        1: reversed_decoy_database,
        2: shuffled_decoy_database,
        3: randomized_decoy_database,
        4: hybrid_decoy_database,
        }
    
    def decoy_database_type(decoy_db_type):
        """
        Decoy Database Types Switcher
        - 1: Reverse Decoy Database
        - 2: Shuffle Decoy Database
        - 3: Random Decoy Database
        - 4: Hybrid Decoy Database
        """
        return decoy_database.decoy_database_switcher.get(decoy_db_type, decoy_database.default)()
    
    def modified_decoy_database(decoy_type):
        modified_decoy_dict = collections.defaultdict(str)
        decoy_list = decoy_database.decoy_database_type(decoy_type)
        for decoy_sequence in decoy_list:
            modified_decoy_dict[decoy_sequence] = decoy_sequence.replace('I','L')
        return modified_decoy_dict
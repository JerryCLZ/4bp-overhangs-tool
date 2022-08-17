from pandas import read_excel
from pandas import DataFrame
from numpy import insert
from itertools import product
from itertools import combinations
from datetime import datetime
from random import randint
from random import sample
from random import uniform
from time import perf_counter
from re import finditer


class Tools:
    
    """This is a class for the tools before designing overhangs
    
    Tools class is the preparations of designing overhangs. It contains a property 
    of DNA codons table, an Ecoli codon usage bias table; general attributes of modules, and files; 
    functions of read_fasta, find_sequence, find_junctions, self_score, overall_score,
    find_position, and seq_bias.
    
    Attributes:
        modules: a list of module names
        files: a string of filenames 
        junc_aa_num: an int of amino acid number
        scores: a table of overhang scores
    """
    
    DNA_codons = {"A":["GCG","GCA","GCC","GCT"],
                  "P":["CCG","CCA","CCC","CCT"],
                  "C":["TGC","TGT"],
                  "Q":["CAG","CAA"],
                  "D":["GAC","GAT"],
                  "E":["GAG","GAA"],
                  "F":["TTC","TTT"],
                  "R":["CGG","CGA","CGC","CGT","AGG","AGA"],
                  "G":["GGG","GGA","GGC","GGT"],
                  "H":["CAC","CAT"],
                  "S":["TCG","TCA","TCC","TCT","AGC","AGT"],
                  "I":["ATA","ATC","ATT"],
                  "T":["ACG","ACA","ACC","ACT"],
                  "K":["AAG","AAA"],
                  "L":["TTG","TTA","CTG","CTA","CTC","CTT"],
                  "V":["GTA","GTC","GTT","GTG"],
                  "W":["TGG"],
                  "Y":["TAC","TAT"],
                  "M":["ATG"],
                  "N":["AAC","AAT"],
                  "Stop":["TAA","TAG","TGA"]
                  }
    
    Ecoli_codon_bias = {'GGG':0.15, 'GGA':0.11, 'GGT':0.34, 'GGC':0.40,
                        'GAG':0.31, 'GAA':0.69, 'GAT':0.63, 'GAC':0.37,
                        'GTG':0.37, 'GTA':0.15, 'GTT':0.26, 'GTC':0.22,
                        'GCG':0.36, 'GCA':0.21, 'GCT':0.16, 'GCC':0.27,
                        'AGG':0.02, 'AGA':0.04, 'CGG':0.10, 'CGA':0.06,
                        'CGT':0.38, 'CGC':0.40, 'AAG':0.23, 'AAA':0.77,
                        'AAT':0.45, 'AAC':0.55, 'ATG':1.00, 'ATA':0.07,
                        'ATT':0.51, 'ATC':0.42, 'ACG':0.27, 'ACA':0.13,
                        'ACT':0.17, 'ACC':0.44, 'TGG':1.00, 'TGT':0.45,
                        'TGC':0.55, 'TAG':0.07, 'TAA':0.64, 'TGA':0.29,
                        'TAT':0.57, 'TAC':0.43, 'TTT':0.57, 'TTC':0.43,
                        'AGT':0.15, 'AGC':0.28, 'TCG':0.15, 'TCA':0.12,
                        'TCT':0.15, 'TCC':0.15, 'CAG':0.65, 'CAA':0.35,
                        'CAT':0.57, 'CAC':0.43, 'TTG':0.13, 'TTA':0.13,
                        'CTG':0.50, 'CTA':0.04, 'CTT':0.10, 'CTC':0.10,
                        'CCG':0.52, 'CCA':0.19, 'CCT':0.16, 'CCC':0.12
                        }
    
    
    def __init__(self, 
                 modules = ('#1: D14_C_chain_A','#2: D14_N_chain_A','#3: D18_C_chain_A'), 
                 files = "example_protein_modules.fasta",
                 junc_aa_num = '2+2',
                 scores_file = '18h_37C',
                 dna_files = "example_DNA_modules.fasta",
                 restriction_file = 'example_restrictions.fasta'
                 ):
        """Initialize Tools class
        
        Input contains list and str
        
        Args:
            modules(list): module names; default is 'D14_chain_A//D14_j1_D79_chain_A//D79_j1_D54_chain_A'
            files(str): module database filenames, use "//"(double forward slash) to separate filenames; default is "test_seq.fasta"
            junc_aa_num(int): the number of amino acid for each junction; default is 6
            scores_file(str): score file name; default is '18h_37C.xlsx'
        """
        self.modules = list(i.split(': ')[1] for i in modules)
        self.files = files
        self.junc_aa_num = junc_aa_num
        self.junc_aa_sum = sum(list(map(int, junc_aa_num.split('+'))))
        self.scores = read_excel(f"{scores_file}.xlsx")
        
        scores_overall = {}
        scores_overhang_axis = self.scores.set_index("Overhang", inplace=False)
        SIV = scores_overhang_axis.index.values
        SCV = scores_overhang_axis.columns.values
        for s1 in SIV:
            for s2 in SCV: #determine one pair of overhangs in the table
                scores_overall[(s1,s2)] = scores_overhang_axis[s1][s2]
        self.scores_overall = scores_overall
        
        scores_self = {}
        scores_number_axis = DataFrame(insert(self.scores.values, 0, values=self.scores.columns, axis=0))
        for x in scores_number_axis.index.values:
            for y in scores_number_axis.index.values:
                if (x==y) & (x!=0):
                    f = scores_number_axis.loc[x,y]
                    s = scores_number_axis.loc[0,x]
                    scores_self[s] = f #grab all diagnal scores from the excel
        self.scores_self = scores_self
        
        max_self_score = max(scores_self.values())
        self.max_self_score = max_self_score
        
        module_protein = {}
        file_list = files.split("\n")
        while "" in file_list:
            file_list.remove("")
        for f in file_list:
            module_lib = open(f) #open files
            dic={}
            for line in module_lib:
                if line.startswith('>'):
                    name=line.replace('>','').strip()
                    dic[name]=''
                else:
                    dic[name]+=line.replace('\n','').strip()  #distribute names to keys and sequences to values
            module_lib.close()
            module_protein = dict(module_protein, **dic)
        self.module_protein = module_protein
        
        module_dna = {}
        dna_file_list = dna_files.split("\n")
        while "" in dna_file_list:
            dna_file_list.remove("")
        for f in dna_file_list:
            module_lib = open(f) #open files
            dic={}
            for line in module_lib:
                if line.startswith('>'):
                    name=line.replace('>','').strip()
                    dic[name]=''
                else:
                    dic[name]+=line.replace('\n','').strip()  #distribute names to keys and sequences to values
            module_lib.close()
            module_dna = dict(module_dna, **dic)
        self.module_dna = module_dna
        
        enzyme_site = {}
        restriction_file_list = restriction_file.split("\n")
        while "" in restriction_file_list:
            restriction_file_list.remove("")
        for f in restriction_file_list:
            module_lib = open(f) #open files
            dic={}
            for line in module_lib:
                if line.startswith('>'):
                    name=line.replace('>','').strip()
                    dic[name]=''
                else:
                    dic[name]+=line.replace('\n','').strip()  #distribute names to keys and sequences to values
            module_lib.close()
            enzyme_site = dict(enzyme_site, **dic)
        self.enzyme_site = enzyme_site
            
    def find_sequence(self, database):
        """This is a sequence searching function
        
        Search the sequences in specified file for the modules
        
        Args:
            None
        
        Returns:
            a dictionary where keys are module names
            and values are sequeneces
        """
        seq = {}
        i = 1
        for mod in self.modules:
            if mod in database:
                seq[str(i)+"_"+mod] = database[mod] #grab the sequences we need from the original dictionary
                i+=1
        return seq

    def find_junctions(self): 
        """This is a junction generating function
        
        Generate the junctions between two sequences, and all possible
        codon combinations for each junction
        
        Args:
            None
        
        Returns:
            a dictionary where keys are junction names and 
            values are their possible codon combinations
        """
        
        def reverse_translate(aa_seq, former_dna, latter_dna, num_base_f, num_base_l, enzyme_site):
            """This is a reverse translation function

            Reverse translate the protein sequence to
            all possible DNA sequences

            Args:
                aa_seq(str): a string of amino acid sequence

            Returns:
                a list containing all possible codon combinations
            """
            enzyme_name_list = []
            site_length_list = []
            overhang_left_list = []
            overhang_right_list = []
            for i in enzyme_site:
                enzyme_name_list.append(i)
                site_length_list.append(len(enzyme_site[i]))
                overhang_left_list.append(former_dna[::-1][num_base_f*3:num_base_f*3+len(enzyme_site[i])-1][::-1])
                overhang_right_list.append(latter_dna[num_base_l*3:num_base_l*3+len(enzyme_site[i])-1])
                        
            seq_codons = []
            for x in aa_seq:
                seq_codons.append(self.DNA_codons[x]) #grab all codons for each amino acid
            result = []
            for r in product(*seq_codons):
                r_joined = "".join(r)
                enzyme_exclude_count = 0
                for i in range(len(enzyme_site)):
                    site_length = site_length_list[i]
                    overhang_left = overhang_left_list[i]
                    overhang_right = overhang_right_list[i]
                    if (enzyme_site[enzyme_name_list[i]] not in overhang_left+r_joined+overhang_right):
                        enzyme_exclude_count += 1
                if enzyme_exclude_count == len(enzyme_site):
                    result.append(r_joined) #generate all possible combinations of these codons
            return result #list
        
        enzyme_site = self.enzyme_site
        num_base_f = int(self.junc_aa_num.split('+')[0])
        num_base_l = int(self.junc_aa_num.split('+')[1])
        protein_seq = self.find_sequence(self.module_protein)
        dna_seq = self.find_sequence(self.module_dna)
        junc_seq = {}
        names = tuple(protein_seq.keys())
        length = len(protein_seq.keys())
        for i in range(length-1):
            former_seq = names[i]
            latter_seq = names[i+1]
            junc_aa = protein_seq[former_seq][::-1][0:num_base_f][::-1] + protein_seq[latter_seq][0:num_base_l]
            junc_seq[str(i+1)+"_"+junc_aa] = reverse_translate(junc_aa, dna_seq[names[i]], dna_seq[names[i+1]], num_base_f, num_base_l, enzyme_site) 
            # find all DNA sequences for each junction
        return junc_seq #dict
    
    def self_score(self, overhang_set):
        """This is a self scoring function

        Based on the combination frequency to score the quality
        of each overhang

        Args:
            overhang_set(list): the overhang(s)

        Returns:
            a list containing the score of each overhang
        """
        score = []
        for o in overhang_set:
            for m in self.scores_self:
                if o[-4:] == m:
                    s = self.scores_self[m]
                    score.append(s/self.max_self_score)
        return score #list, return the corresponding scores
    
    def overall_score(self, overhang_set):
        """This is a overall scoring function
        
        Based on the combination frequency to score the quality
        of the set of overhangs
        
        Args:
            overhang_set(list): the overhangs (no less than 2)
            
        Returns:
            the score which sums up all the scores between each overhang
        """ 
                        
        def reverse_complementary(seq):
            """This is a reverse and complementary sequence generating function

            Based on the principle of complementary base pairing, generate 
            the complementary sequence of the given DNA sequence

            Args:
                seq(str): a string of DNA sequence

            Returns:
                a string of the complementary sequence or the error message
            """
            base_paring = {"A":"T","T":"A","C":"G","G":"C"}
            result = ""
            for n in seq:
                for b in base_paring:
                    if n==b:
                        result += base_paring[b]
            return result[::-1]
        
        maxmax = self.max_self_score**2
        scores = self.scores_overall
        score = 0
        overhang_from_set = list(combinations(overhang_set, 2)) #determine one pair of overhangs in the set
        cbn_score = {} #record the score of each combination
        for os in overhang_from_set:
            if os[0][-4:] == os[1][-4:]:
                score = maxmax
                break
            os1 = (os[0][-4:],reverse_complementary(os[1][-4:]))
            os2 = (reverse_complementary(os[0][-4:]),os[1][-4:])
            score += scores[os1]+scores[os2]
            cbn_score[os] = scores[os1]+scores[os2]
            #if two pairs are same, or not fully complementary,
            #they are the possible mismatiched overhangs and the score should be added
        for o in overhang_set:
            o_score = 0
            for c in cbn_score:
                if o in c:
                    o_score += cbn_score[c]
            if o_score > 0.5*(self.self_score([o])[0]):
                score = maxmax
                break         
        # to prevent there is any overhang which has high mismatch score
        return round(score*2/(len(overhang_set)**2-len(overhang_set)), 4)
    
    def find_position(self, overhang_set, Junctions):
        """This is a overhang postition searching function
        
        Find the all postions of the first base of the overhangs
        in their own junctions
        
        Args:
            overhang_set(list):the overhang(s)
            Junctions(dict):a dictionary contains junction sequences and their DNA sequences
            
        Returns:
            a dictionary where the keys are overhangs and 
            values are their positions
        """
        if overhang_set == []:
            overhang_position = {'0_0':[0]}
        else:
            emptyp = []
            for i in range(len(overhang_set)):
                emptyp.append([]) #generate same length list, to store the position number
            overhang_position = dict(zip(overhang_set,emptyp)) #combine the overhang list and empty position list together
            for junc_name in Junctions:
                junc_num = tuple(Junctions.keys()).index(junc_name)
                overhang = overhang_set[junc_num]
                for s in range(len(Junctions[junc_name])):
                    if overhang[-4:] in Junctions[junc_name][s]:
                        position = []
                        for p in finditer(overhang[-4:], Junctions[junc_name][s]):
                            position.append(p.start())
                        if position not in overhang_position[overhang]:
                            overhang_position[overhang].append(position)
        return overhang_position
    
    def most_least_common_codons(self):
        """This is a most commonly used codons generating function
        
        Based on Ecoli codon usage bias table, to generate the most commonly used codon for each amino acid
        
        Args:
            None
            
        Returns:
            a dictionary where keys are amino acids and values are the frequency of the most common one
        """
        Ecoli_most_common_codons = {}
        for aa in self.DNA_codons:
            max_bias = 0
            for codon in self.DNA_codons[aa]:
                bias = self.Ecoli_codon_bias[codon]
                if bias > max_bias:
                    max_bias = bias
            Ecoli_most_common_codons[aa] = max_bias
        
        Ecoli_least_common_codons = {}
        for aa in self.DNA_codons:
            min_bias = 1
            for codon in self.DNA_codons[aa]:
                bias = self.Ecoli_codon_bias[codon]
                if bias < min_bias:
                    min_bias = bias
            Ecoli_least_common_codons[aa] = min_bias
        return Ecoli_most_common_codons, Ecoli_least_common_codons
    
    def seq_bias(self, overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons):
        """This is a function for searching the junction DNA sequence with the highest bias score 
        
        Based on the overhangs and positions, find all possible DNA sequences in each junctions, 
        then return the one with highest bias score
        
        Args:
            overhang_position(dict): a dictionary contains overhangs and their positions
            Junctions(dict):a dictionary contains junction sequences and their DNA sequences
            Ecoli_most_common_codons(dict): a dictionary containing amino acids and their mostly used codons
        
        Returns:
            a dictionary where keys are junction DNA sequences and values are their codon bias percentage
        """
        
        def bias_score(seq, junction, Ecoli_most_common_codons, Ecoli_least_common_codons):
            """This is a codon bias scoring function
            
            Based on the E.coli codon usage bias, add up the frequency of all codons in the sequnece
            
            Args:
                seq(str): a string of DNA sequence translated from amino acid sequence
                junction(str): a string of junction in amino acid sequence
                Ecoli_most_common_codons(dict): a dictionary containing amino acids and their mostly used codons
                
            Returns:
                a float number, the ratio of the sum of bias scores to the sum of max scores
            """
            if len(seq)%3 != 0:
                return "Not divisible by 3"
            else:
                min_score = 0
                for a in junction:
                    min_score += Ecoli_least_common_codons[a]
                max_score = 0
                for a in junction:
                    max_score += Ecoli_most_common_codons[a]
                score = 0
                for i in range(0,len(seq),3):
                    score += self.Ecoli_codon_bias[seq[i:i+3]]
                return round((score-min_score)/(max_score-min_score), 2)
        
        if overhang_position == {'0_0':[0]}:
            best_oseq = {'0_0':0}
        else:
            best_oseq = {}
            sum_score = 0
            for junc_name in Junctions:
                oseq_score = {}
                num = tuple(Junctions.keys()).index(junc_name)
                overhang = tuple(overhang_position.keys())[int(num)]
                position = overhang_position[overhang]
                for seq in Junctions[junc_name]:
                    if overhang[-4:] in seq:
                        if [p.start() for p in finditer(overhang[-4:],seq)] in position:
                            score = bias_score(seq, junc_name.split("_")[1], Ecoli_most_common_codons, Ecoli_least_common_codons)
                            oseq_score[seq] = score
                oseq_score_sorted = sorted(oseq_score.items(), key=lambda x:x[1], reverse=True)
                best_oseq[str(num+1)+"_"+oseq_score_sorted[0][0]] = oseq_score_sorted[0][1]
        return best_oseq

    def find_primer(self, best_overhang_position, max_seq_bias):
        """
        best_overhang_position(dict)
        max_seq_bias(dict)
        """
        def Dot_Seq(DNA_sequence, junc_sequence, former_len, offset, direction):
            junc_sequence = junc_sequence.split('_')[1]
            if direction == 'front':
                dot_seq = ''
                if DNA_sequence[offset+13] in ['C','G']:
                    dot_seq = DNA_sequence[offset:offset+14]
                if DNA_sequence[offset+14] in ['C','G']:
                    dot_seq = DNA_sequence[offset:offset+15]
                if DNA_sequence[offset+15] in ['C','G']:
                    dot_seq = DNA_sequence[offset:offset+16]
                if dot_seq == '':
                    dot_seq = DNA_sequence[offset:offset+15]+'C'        
            if direction == 'back':
                dot_seq = ''
                if DNA_sequence[offset-14] in ['C','G']:
                    dot_seq = DNA_sequence[::-1][former_len:14-offset][::-1]+junc_sequence[:(former_len+offset)]
                if DNA_sequence[offset-15] in ['C','G']:
                    dot_seq = DNA_sequence[::-1][former_len:15-offset][::-1]+junc_sequence[:(former_len+offset)]
                if DNA_sequence[offset-16] in ['C','G']:
                    dot_seq = DNA_sequence[::-1][former_len:16-offset][::-1]+junc_sequence[:(former_len+offset)]
                if dot_seq == '':
                    dot_seq = 'C'+DNA_sequence[::-1][former_len:15-offset][::-1]+junc_sequence[:(former_len+offset)]
            return dot_seq
        
        if (best_overhang_position == {'0_0':[0]}) or (max_seq_bias == {'0_0':0}):
            Primer_list = {'0_0':0}
        else:
            Primer_list = {}
            dna_seq = self.find_sequence(self.module_dna)
            for i in dna_seq:
                if int(i.split('_')[0]) == 1:
                    position = list(best_overhang_position.values())[int(i.split('_')[0])-1][0]
                    former_len = list(map(int, self.junc_aa_num.split('+')))[0]*3
                    Primer_list[i+'_US'] = 'TCAGCATATG' + Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, 0, 'front')       
                    if position < former_len-3:
                        extra_seq = ''
                        offset = -(former_len-position-4)
                        dot_seq = Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, offset, 'back')
                    if position >= former_len-3:
                        extra_seq = list(max_seq_bias.keys())[int(i.split('_')[0])-1].split('_')[1][former_len:position+4]
                        dot_seq = Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, 0, 'back')
                    Primer_list[i+'_DS'] = dot_seq + extra_seq + 'TGAGACCCTCGAGTAA'
                elif int(i.split('_')[0]) == len(dna_seq):
                    former_len = list(map(int, self.junc_aa_num.split('+')))[0]*3
                    Primer_list[i+'_US'] = 'TCAGCATATGAGGTCTCC' + Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[len(max_seq_bias)-1], 0, 0, 'front')
                    Primer_list[i+'_DS'] = Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[len(max_seq_bias)-1], 0, 0, 'back') + 'GGTTGGCTCGAGATAG'
                    # list(max_seq_bias.keys())[len(max_seq_bias)-1] is useless actually, just to keep the input format.
                else:
                    position = list(best_overhang_position.values())[int(i.split('_')[0])-1][0]
                    former_len = list(map(int, self.junc_aa_num.split('+')))[0]*3
                    Primer_list[i+'_US'] = 'TCAGCATATGAGGTCTCC' + Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, 0, 'front')
                    if position < former_len-3:
                        extra_seq = ''
                        offset = -(former_len-position-4)
                        dot_seq = Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, offset, 'back')
                    if position >= former_len-3:
                        extra_seq = list(max_seq_bias.keys())[int(i.split('_')[0])-1].split('_')[1][former_len:position+4]
                        dot_seq = Dot_Seq(dna_seq[list(dna_seq.keys())[int(i.split('_')[0])-1]], list(max_seq_bias.keys())[int(i.split('_')[0])-1], former_len, 0, 'back')
                    Primer_list[i+'_DS'] = dot_seq + extra_seq + 'TGAGACCCTCGAGTAA'
        return Primer_list       

class Operators:
    """This is a class for genetic algorithm

    Operators class is the set of arguments and functions of genetic algorithm. 
    It contains general attributes of init_num, mutation_rate, crossover_rate, 
    generation_num, and pre; and functions of start, select, mutate, and crossover.

    Attributes:
        init_num: an int number of initial individuals
        mutation_rate: a float number of mutation rate
        crossover_rate: a float number of crossover rate
        pre: a object of a class; 
    """
    
    def __init__(self, init_num=20, mutation_rate=0.02, crossover_rate=0.2, pre=Tools()):

        """Initialize Operators class
        
        details
        
        Args:
            init_num(int): initial number of individuals, must be divisible by 10; default is 50
            mutation_rate(float): mutation rate; default is 0.01
            crossover_rate(float): crossover rate; default is 0.2
            pre(object): object name; default is the objcet of class Tools
        """
        self.init_num = init_num
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.pre = pre
    
    def start(self, Junctions, init_num):
        """This is a initial individuals generating function
        
        Generate a certain number(init_num) of initial individuals randomly.
        Each individual contains an overhang set and its overall score.
        
        Args:
            Junctions(dict): a dictionary of the junction names and their sequences
            init_num(number): an int number to determine the number of initial individuals
            
        Returns:
            a dictionary where the keys are overhang sets and 
            values are their overall score
        """
        individuals = {}
        for num in range(init_num):
            overhang_set = []
            for junc_name in Junctions:
                junc_num = tuple(Junctions.keys()).index(junc_name)
                p = randint(0 , self.pre.junc_aa_sum*3-4)
                seq = sample(list(Junctions[junc_name]),1)
                overhang = seq[0][p:p+4]
                overhang_set.append(str(junc_num+1)+"_"+overhang)
            ov_score = self.pre.overall_score(overhang_set)
            individuals[" ".join(overhang_set)] = ov_score
        return individuals  
        
    def select(self, indiv_set, self_score_lower_limit): #what to do if most or all individuals contain the <500 overhang?
        """This is a selection function
        
        Filter out the individuals which contains the overhang below a certain self score; then filter
        the top half individuals based on overall score
        
        Args:
            indiv_set(dict): a dictionary where keys are overhang sets and values are their overall scores
            self_score_lower_limit(int): a number that the lower limit of self score is
        Returns:
            a dictionary where keys are overhang sets and values are their overall scores
        """
        num = len(indiv_set)*0.5
        indiv_set_HQ = dict(tuple(indiv_set.items()))
        for i in indiv_set:
            o = i.split(" ")
            scores = self.pre.self_score(o)
            for s in scores:
                if s < self_score_lower_limit:
                    del indiv_set_HQ[i]
                    break
        indiv_set_sorted = sorted(indiv_set_HQ.items(), key=lambda x:x[1], reverse=False)
        tops = dict(indiv_set_sorted[0:int(num)])
        return tops
        
    def mutate(self, indiv_set, Junctions):
        """This is a mutation function
        
        Each overhang in each individuals has a certain chance to be replaced by a random
        overhang from the same junction pool
        
        Args:
            indiv_set(list): a list of individuals
            Junctions(dict): a dictionary containing all DNA sequences in each junction
        
        Returns:
            a list of the mutated individuals
        """
        for n in range(len(indiv_set)):
            overhang_set = indiv_set[n].split()
            for i in range(len(overhang_set)):
                x = uniform(0,1)
                if x < self.mutation_rate:
                    junc_name = tuple(Junctions.keys())[i]
                    p = randint(0 , self.pre.junc_aa_sum*3-4)
                    seq = sample(list(Junctions[junc_name]),1)
                    overhang = seq[0][p:p+4]
                    overhang_set[i] = str(i+1)+"_"+overhang
            indiv_set[n] = " ".join(overhang_set)
        return indiv_set

    def crossover(self, indiv_set):
        """This is a crossover function
        
        Grouping all individuals in random pairs. for each junction, the paired overhang sets
        has a certain chance to exchange their overhangs.
        
        Args:
            indiv_set(list): a list of individuals
        
        Returns:
            a list of crossovered individuals
        """
        if indiv_set == []:
            return indiv_set
        else:
            indiv_set_listed = []
            for i in indiv_set:
                indiv_set_listed.append(i.split())
            junc_num = len(indiv_set[0].split())

            crossed_set = []
            cross_group = []
            for g in range(int(len(indiv_set)/2)):
                group_num = sample(range(len(indiv_set_listed)), 2)
                group_num = sorted(group_num,reverse=False)
                group = []
                for a in list(group_num):
                    group.append(indiv_set_listed[a])
                cross_group.append(group)
                indiv_set_listed_dep = []
                for y in range(len(indiv_set_listed)):
                    if y not in group_num:
                        indiv_set_listed_dep.append(indiv_set_listed[y])
                indiv_set_listed = indiv_set_listed_dep
            for n in range(len(cross_group)):
                grps = cross_group[n]
                for j in range(junc_num):
                    x = uniform(0,1)
                    if x < self.crossover_rate:
                        dep = grps[0][j]
                        grps[0][j] = grps[1][j]
                        grps[1][j] = dep
                crossed_set.append(" ".join(grps[0]))
                crossed_set.append(" ".join(grps[1]))
            return crossed_set

# Monte carlo
def Montecarlo(files, DNA_files, Restriction_file, modules, times_uppser_limit, times_lower_limit, self_score_lower_limit,
               aa_num, score_table, loadingwin, loadingbar, Ecoli_codon_bias):
    
    time_start = perf_counter() #timing
    monte_pre = Tools(modules,files,aa_num,score_table, DNA_files, Restriction_file) #tools object
    Junctions = monte_pre.find_junctions() #generate all possible sequences in the jucntions
    Ecoli_most_common_codons, Ecoli_least_common_codons = monte_pre.most_least_common_codons()

    min_ov_score = 9999
    times = 0
    max_seq_bias = {"":0}
    best_overhang_position = []

    while (min_ov_score > 0) or (sum(max_seq_bias.values())/len(Junctions) < Ecoli_codon_bias) or (times < times_lower_limit):
        #stop condition is that overall score is 0;
        overhang_set = []
        for junc_name in Junctions:
            junc_num = tuple(Junctions.keys()).index(junc_name)
            overhang_score = self_score_lower_limit
            while overhang_score <= self_score_lower_limit:
                p = randint(0 , monte_pre.junc_aa_sum*3-4)
                seq = sample(list(Junctions[junc_name]),1)
                overhang = seq[0][p:p+4]
                overhang_score = monte_pre.self_score([overhang])[0]
            overhang_set.append(str(junc_num+1)+"_"+overhang)
            # randomly generate the overhang whose self score is above self_score_lower_limit
        ov_score = monte_pre.overall_score(overhang_set) #calculate overall score
        if ov_score < min_ov_score:
            overhang_position = monte_pre.find_position(overhang_set, Junctions)
            overhang_seq_bias= monte_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons)
            min_ov_score = ov_score
            max_seq_bias = overhang_seq_bias
            best_overhang_position = overhang_position
        if Ecoli_codon_bias == 1:
            if ov_score == min_ov_score:
                overhang_position = monte_pre.find_position(overhang_set, Junctions)
                overhang_seq_bias= monte_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons)
                if sum(max_seq_bias.values()) <= sum(overhang_seq_bias.values()):
                    min_ov_score = ov_score
                    max_seq_bias = overhang_seq_bias
                    best_overhang_position = overhang_position
        # record the overhang set with minimal overall score and highest bias score
        loadingbar.step(1)
        loadingwin.update()
        times += 1
        if times >= times_uppser_limit:
            break # if exceed max cycle times, then stop
    for o,s in zip(best_overhang_position, max_seq_bias):
        best_overhang_position[o] = [p.start() for p in finditer(o[-4:], s.split("_")[1])]
    time_end = perf_counter() #timing

    OFF_ON=["OFF","ON"]
    aa_num_former = list(map(int, aa_num.split('+')))[0]
    aa_sum = sum(list(map(int, aa_num.split('+'))))
    juncs = "  ".join(Junctions.keys())
    filetime = datetime.now().strftime("%c").replace(":","-").replace(" ","_")
    primer_list = monte_pre.find_primer(best_overhang_position, max_seq_bias)
    with open(f"Monte Carlo Running Report - {filetime}.txt", "w") as f:
        f.write(f"This is the result from Monte Carlo algorithm\n\n")
        f.write(f"Mostly used codons: {OFF_ON[Ecoli_codon_bias]}\nResidue number: {aa_num}\nLigation time: {score_table.split('_')[0]}\nTempreture: {score_table.split('_')[1]}\n")
        f.write(f"Overhang quality limit: {self_score_lower_limit}\nMaximum cycles: {times_uppser_limit}\nMinimum cycles: {times_lower_limit}\n\n")
        f.write(f"Protein files:\n{files}\n\nDNA files:\n{DNA_files}\n\n")
        f.write(f"Restriction sites avoided: ")
        for i in monte_pre.enzyme_site:
            f.write(f"{i}, ")
        f.write(f"\n{Restriction_file}\n\n")
        f.write(f"Modules:\n")
        for i in modules:
            f.write(f"{i}  ")
        f.write(f"\n\nJunctions:\n{juncs}\n\n")
        f.write(f"totally {times} cycles\nRunning time: {time_end-time_start} seconds\n\n")
        f.write(f"the overall mismatch score is {min_ov_score}\n")
        f.write(f"No.  Overhang  former3'end(1-{aa_num_former*3})  latter5'end({aa_num_former*3+1}-{aa_sum*3})  Position(1-{aa_sum*3})  Codon_usage_score\n")
        for a,b in zip(best_overhang_position, max_seq_bias):
            col_No = int(a.split("_")[0])
            col_Overhang = a.split("_")[1]
            col_DNA_Sequence = b.split("_")[1]
            col_Position = ", ".join(list(map(lambda x:str(x+1),best_overhang_position[a])))
            col_Codon_bias_score = max_seq_bias[b]
            f.write(f"{col_No}    {col_Overhang}    {col_DNA_Sequence[0:aa_num_former*3]}    {col_DNA_Sequence[aa_num_former*3:]}    {col_Position}    {col_Codon_bias_score}\n")
        f.write(f"\nPrimer details:\n")
        for i in primer_list:
            f.write(f"{i}: {primer_list[i]}\n")   
    loadingbar['value'] = times_uppser_limit+1
            
# Greedy algorithm
def Greedy(files, DNA_files, Restriction_file, modules, self_score_lower_limit, aa_num, score_table, loadingwin, loadingbar, Ecoli_codon_bias):
    
    time_start = perf_counter() # timing
    aa_sum = sum(list(map(int, aa_num.split('+'))))
    greedy_pre = Tools(modules,files,aa_num,score_table, DNA_files, Restriction_file) # tools object
    Junctions = greedy_pre.find_junctions()# generate all sequences in the Junctions
    Ecoli_most_common_codons, Ecoli_least_common_codons = greedy_pre.most_least_common_codons()
    overhang_set = []

    for junc_name in Junctions:
        junc_num = tuple(Junctions.keys()).index(junc_name)
        indiv_scores = {}
        for s in range(len(Junctions[junc_name])):
            for p in range(0,aa_sum*3-4):
                overhang = Junctions[junc_name][s][p:p+4]
                score = greedy_pre.self_score([overhang])[0]
                if overhang not in indiv_scores:
                    if score > self_score_lower_limit:
                        indiv_scores[overhang] = score
        indiv_scores = dict(sorted(indiv_scores.items(), key=lambda x:x[1], reverse=True))
        #for each junction, generate and sort all possible overhangs and their self score
        if len(overhang_set) == 0:
            highest_bias = 0
            highest_bias_position = 0
            bias = 0
            position = 0
            overhang = tuple(indiv_scores.keys())[position] #if codon bias=0, then output the first one; if 1, then go into loop.
            if Ecoli_codon_bias == 1:
                while bias < Ecoli_codon_bias:
                    if position >= len(indiv_scores):
                        overhang = tuple(indiv_scores.keys())[highest_bias_position]
                        break
                    overhang = tuple(indiv_scores.keys())[position]
                    junc_name = tuple(Junctions.keys())[len(overhang_set)]
                    junc_seq = tuple(Junctions.values())[len(overhang_set)]
                    overhang_p = greedy_pre.find_position([overhang], dict(zip([junc_name], [junc_seq])))
                    bias = tuple(greedy_pre.seq_bias(overhang_p, dict(zip([junc_name], [junc_seq])), Ecoli_most_common_codons, Ecoli_least_common_codons).values())[0]
                    if bias > highest_bias:
                        highest_bias = bias
                        highest_bias_position = position
                    position += 1
            overhang_set.append("1_"+overhang)
        # the first one will be recorded directly, after attaining the criterias
        else:
            ov_score = 1
            lowest_ov = max(greedy_pre.scores_overall.values())
            lowest_ov_position = 0
            bias = 0
            position = 0
            while (ov_score > 0) or (bias < Ecoli_codon_bias):
                if position >= len(indiv_scores):
                    overhang = tuple(indiv_scores.keys())[lowest_ov_position]
                    break
                overhang = tuple(indiv_scores.keys())[position]
                ov_score = greedy_pre.overall_score(overhang_set+[overhang])
                junc_name = tuple(Junctions.keys())[len(overhang_set)]
                junc_seq = tuple(Junctions.values())[len(overhang_set)]
                overhang_p = greedy_pre.find_position([overhang], dict(zip([junc_name], [junc_seq])))
                bias = tuple(greedy_pre.seq_bias(overhang_p, dict(zip([junc_name], [junc_seq])), Ecoli_most_common_codons, Ecoli_least_common_codons).values())[0]
                if ov_score < lowest_ov:
                    lowest_ov = ov_score
                    lowest_ov_position = position
                if Ecoli_codon_bias == 1:
                    if ov_score == lowest_ov:
                        if bias > highest_bias:
                            highest_bias = bias
                            lowest_ov_position = position
                position += 1
            overhang_set.append(str(junc_num+1)+"_"+overhang)
        loadingbar.step(1)
        loadingwin.update()
        # the followings should be checked for overall score, if not attained, use the next one in this junction; if runs out, return "None"
    ov_score = greedy_pre.overall_score(overhang_set)
    if ov_score > 9999:
        ov_score = 9999
    overhang_position = greedy_pre.find_position(overhang_set, Junctions) #find the position
    seq_bias = greedy_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons) #find junction sequences
    for o,s in zip(overhang_position, seq_bias):
        overhang_position[o] = [p.start() for p in finditer(o[-4:], s.split("_")[1])]
    time_end = perf_counter()#timing

    OFF_ON=["OFF","ON"]
    aa_num_former = list(map(int, aa_num.split('+')))[0]
    juncs = "  ".join(Junctions.keys())
    filetime = datetime.now().strftime("%c").replace(":","-").replace(" ","_")
    primer_list = greedy_pre.find_primer(overhang_position, seq_bias)
    with open(f"Greedy Running Report - {filetime}.txt", "w") as f:        
        f.write(f"This is the result from Greedy algorithm\n\n")
        f.write(f"Mostly used codons: {OFF_ON[Ecoli_codon_bias]}\nResidue number: {aa_num}\nLigation time: {score_table.split('_')[0]}\nTempreture: {score_table.split('_')[1]}\n")
        f.write(f"Overhang quality limit: {self_score_lower_limit}\n\n")
        f.write(f"Protein files:\n{files}\n\nDNA files:\n{DNA_files}\n\n")
        f.write(f"Restriction sites avoided: ")
        for i in greedy_pre.enzyme_site:
            f.write(f"{i}, ")
        f.write(f"\n{Restriction_file}\n\n")
        f.write(f"Modules:\n")
        for i in modules:
            f.write(f"{i}  ")
        f.write(f"\n\nJunctions:\n{juncs}\n\n")
        f.write(f"Running time: {time_end-time_start} seconds\n\n")        
        f.write(f"the overall mismatch score is {ov_score}\n")
        f.write(f"No.  Overhang  former3'end(1-{aa_num_former*3})  latter5'end({aa_num_former*3+1}-{aa_sum*3})  Position(1-{aa_sum*3})  Codon_usage_score\n")
        for a,b in zip(overhang_position, seq_bias):
            col_No = int(a.split("_")[0])
            col_Overhang = a.split("_")[1]
            col_DNA_Sequence = b.split("_")[1]
            col_Position = ", ".join(list(map(lambda x:str(x+1),overhang_position[a])))
            col_Codon_bias_score = seq_bias[b]
            f.write(f"{col_No}    {col_Overhang}    {col_DNA_Sequence[0:aa_num_former*3]}    {col_DNA_Sequence[aa_num_former*3:]}    {col_Position}    {col_Codon_bias_score}\n")
        f.write(f"\nPrimer details:\n")
        for i in primer_list:
            f.write(f"{i}: {primer_list[i]}\n") 
    loadingbar['value'] = len(modules)+1

# Genetic algorithm
def Genetic(files, DNA_files, Restriction_file, modules, init_num, mutation_rate, crossover_rate, max_generation_num, min_generation_num,
            self_score_lower_limit, aa_num, score_table, loadingwin, loadingbar, Ecoli_codon_bias):

    time_start = perf_counter() #timing
    genetic_pre = Tools(modules,files,aa_num,score_table, DNA_files, Restriction_file) #tools object
    operators = Operators(init_num=init_num, 
                          mutation_rate=mutation_rate, 
                          crossover_rate=crossover_rate,
                          pre = genetic_pre) #operators object
    Junctions = genetic_pre.find_junctions() #generate all possible DNA sequences in the junctions
    Ecoli_most_common_codons, Ecoli_least_common_codons = genetic_pre.most_least_common_codons()
    best_ov_score = 9999
    best_average_bias = {"":0}
    best_indiv = ''

    initials = operators.start(Junctions, operators.init_num) #generate the initial individuals

    G = 0
    while (best_ov_score > 0) or (sum(best_average_bias.values())/len(Junctions) < Ecoli_codon_bias) or (G < min_generation_num):
        indiv_set_selected = operators.select(initials, self_score_lower_limit) #filter the initials by overall score and self score
        indiv_set_mutated = operators.mutate(list(indiv_set_selected.keys()), Junctions) #mutate them
        indiv_set_crossed = operators.crossover(indiv_set_mutated) #crossover them
        initials2 = {}
        for overhang_set in indiv_set_crossed:
            overhang_set = overhang_set.split()
            ov_score = genetic_pre.overall_score(overhang_set)
            initials2[" ".join(overhang_set)] = ov_score
        initials2 = dict(list(initials2.items())+list(indiv_set_selected.items())) #put the operated individuals and originals together
        indiv_set_selected2 = operators.select(initials2, self_score_lower_limit) #filter the operateds and olds together
        for indiv, score in indiv_set_selected2.items():
            if score < best_ov_score:
                overhang_position = genetic_pre.find_position(indiv.split(), Junctions)
                bias = genetic_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons)
                best_ov_score = score
                best_average_bias = bias
                best_indiv = indiv  # record the best individual
            if Ecoli_codon_bias == 1:
                if score == best_ov_score:
                    overhang_position = genetic_pre.find_position(indiv.split(), Junctions)
                    bias = genetic_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons)
                    if sum(bias.values()) >= sum(best_average_bias.values()): 
                        best_ov_score = score
                        best_average_bias = bias
                        best_indiv = indiv  # record the best individual
        new_indiv_set = operators.start(Junctions, (operators.init_num - len(indiv_set_selected2))) #add new random individuals to the initial number
        indiv_set_next_generation = dict(list(indiv_set_selected2.items())+list(new_indiv_set.items())) # put old ones and new ones together
        initials = {}
        for overhang_set in indiv_set_next_generation:
            overhang_set = overhang_set.split()
            ov_score = genetic_pre.overall_score(overhang_set)
            initials[" ".join(overhang_set)] = ov_score
        #regard the new generation as new initials
        loadingbar.step(1)
        loadingwin.update()
        G += 1
        if G >= max_generation_num:
            break

    overhang_position = genetic_pre.find_position(best_indiv.split(), Junctions)
    seq_bias = genetic_pre.seq_bias(overhang_position, Junctions, Ecoli_most_common_codons, Ecoli_least_common_codons)
    for o,s in zip(overhang_position, seq_bias):
        overhang_position[o] = [p.start() for p in finditer(o[-4:], s.split("_")[1])]
    time_end = perf_counter() #timing

    OFF_ON=["OFF","ON"]
    aa_num_former = list(map(int, aa_num.split('+')))[0]
    aa_sum = sum(list(map(int, aa_num.split('+'))))
    juncs = "  ".join(Junctions.keys())
    filetime = datetime.now().strftime("%c").replace(":","-").replace(" ","_")
    primer_list = genetic_pre.find_primer(overhang_position, seq_bias)
    with open(f"Genetic Running Report - {filetime}.txt", "w") as f:
        f.write(f"This is the result from Genetic algorithm\n\n")
        f.write(f"Mostly used codons: {OFF_ON[Ecoli_codon_bias]}\nResidue number: {aa_num}\nLigation time: {score_table.split('_')[0]}\nTempreture: {score_table.split('_')[1]}\n")
        f.write(f"Overhang quality limit: {self_score_lower_limit}\nMaximum generations: {max_generation_num}\nMinimum generations: {min_generation_num}\nInitial number: {init_num}\nMutation rate: {mutation_rate}\nCrossover rate: {crossover_rate}\n\n")
        f.write(f"Protein files:\n{files}\n\nDNA files:\n{DNA_files}\n\n")
        f.write(f"Restriction sites avoided: ")
        for i in genetic_pre.enzyme_site:
            f.write(f"{i}, ")
        f.write(f"\n{Restriction_file}\n\n")
        f.write(f"Modules:\n")
        for i in modules:
            f.write(f"{i}  ")
        f.write(f"\n\nJunctions:\n{juncs}\n\n")
        f.write(f"totally {G} generations\nRunning time: {time_end-time_start} seconds\n\n")
        f.write(f"the overall mismatch score is {best_ov_score}\n")
        f.write(f"No.  Overhang  former3'end(1-{aa_num_former*3})  latter5'end({aa_num_former*3+1}-{aa_sum*3})  Position(1-{aa_sum*3-3})  Codon_usage_score\n")
        for a,b in zip(overhang_position, seq_bias):
            col_No = int(a.split("_")[0])
            col_Overhang = a.split("_")[1]
            col_DNA_Sequence = b.split("_")[1]
            col_Position = ", ".join(list(map(lambda x:str(x+1),overhang_position[a])))
            col_Codon_bias_score = seq_bias[b]
            f.write(f"{col_No}    {col_Overhang}    {col_DNA_Sequence[0:aa_num_former*3]}    {col_DNA_Sequence[aa_num_former*3:]}    {col_Position}    {col_Codon_bias_score}\n")
        f.write(f"\nPrimer details:\n")
        for i in primer_list:
            f.write(f"{i}: {primer_list[i]}\n") 
    loadingbar['value'] = max_generation_num+1
    
    
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
from windnd import hook_dropfiles
from os.path import isfile

aa_num_list = ['1+1','2+2','3+3','0+1','0+2','0+3','1+0','2+0','3+0']

def rightside():
    if var1.get()==0:
        lbIN.place_forget()
        textIN.place_forget()
        lbMR.place_forget()
        textMR.place_forget()
        lbCR.place_forget()
        textCR.place_forget()
        lbGUL.place_forget()
        textGUL.place_forget()
        lbGLL.place_forget()
        textGLL.place_forget()
        lbCUL.place(x=660,y=120,height=30,width=200)
        textCUL.place(x=860,y=120,height=30,width=80)
        lbCLL.place(x=660,y=160,height=30,width=200)
        textCLL.place(x=860,y=160,height=30,width=80)
        biascheck.select()
    if var1.get()==1:
        lbIN.place_forget()
        textIN.place_forget()
        lbMR.place_forget()
        textMR.place_forget()
        lbCR.place_forget()
        textCR.place_forget()
        lbGUL.place_forget()
        textGUL.place_forget()
        lbGLL.place_forget()
        textGLL.place_forget()
        lbCUL.place_forget()
        textCUL.place_forget()
        lbCLL.place_forget()
        textCLL.place_forget()
        biascheck.select()
    if var1.get()==2:
        lbCUL.place_forget()
        textCUL.place_forget()
        lbCLL.place_forget()
        textCLL.place_forget()
        lbIN.place(x=660,y=200,height=30,width=200)
        textIN.place(x=860,y=200,height=30,width=80)
        lbMR.place(x=660,y=240,height=30,width=200)
        textMR.place(x=860,y=240,height=30,width=80)
        lbCR.place(x=660,y=280,height=30,width=200)
        textCR.place(x=860,y=280,height=30,width=80)
        lbGUL.place(x=660,y=120,height=30,width=200)
        textGUL.place(x=860,y=120,height=30,width=80)
        lbGLL.place(x=660,y=160,height=30,width=200)
        textGLL.place(x=860,y=160,height=30,width=80)
        biascheck.select()
    
def dragged_DBfiles(files):
    msg = '\n'.join((item.decode('gbk') for item in files))+'\n'
    DBtext.insert(tk.END,msg)

def dragged_LRfiles(files):
    msg = '\n'.join((item.decode('gbk') for item in files))+'\n'
    textLR.insert(tk.END,msg)

def dragged_DNAfiles(files):
    msg = '\n'.join((item.decode('gbk') for item in files))+'\n'
    DNAtext.insert(tk.END,msg)
    
def dragged_RSfiles(files):
    msg = '\n'.join((item.decode('gbk') for item in files))+'\n'
    RStext.insert(tk.END,msg)
    
def load_files():
    MDbox1.delete(0,tk.END)
    MDbox2.delete(0,tk.END)
    files = DBtext.get('1.0',tk.END+"-2c")
    if files != "":
        file_list = files.split("\n")
        while "" in file_list:
            file_list.remove("")
        module_list = []
        errors = ""
        for file in file_list:
            if isfile(file) == True:
                module_lib = open(file) #open files
                for line in module_lib:
                    if line.startswith('>'):
                        name=line.replace('>','').strip()
                        module_list.append(name)
                module_lib.close()
            else:
                errors += f"{file}\n"
        for m in module_list:
            MDbox1.insert(tk.END, m)
        if errors != "":
            tk.messagebox.showerror('Error',f'{errors}do not exist')

def add_module():
    module = MDbox1.get(MDbox1.curselection())
    length = MDbox2.size()+1
    MDbox2.insert(tk.END, f"#{length}: {module}")
    
def delete_module():
    MDbox2.delete(tk.END)

def clear_module():
    MDbox2.delete(0,tk.END)
    
def search_module():
    target = SearchE.get()
    if target == "":
        MDbox1.delete(0,tk.END)
        load_files()
    else:
        module_list = []
        for module in MDbox1.get(0,tk.END):
            if target in module:
                module_list.append(module)
        if module_list != []:
            MDbox1.delete(0,tk.END)
            for m in module_list:
                MDbox1.insert(tk.END, m)
    
def clearall():
    textSSLL.delete('1.0',tk.END)
    textSSLL.insert(tk.END, '0.1')
    textCUL.delete('1.0',tk.END)
    textCUL.insert(tk.END, '10000')
    textCLL.delete('1.0',tk.END)
    textCLL.insert(tk.END, '100')
    textIN.delete('1.0',tk.END)
    textIN.insert(tk.END, '10')
    textMR.delete('1.0',tk.END)
    textMR.insert(tk.END, '0.1')
    textCR.delete('1.0',tk.END)
    textCR.insert(tk.END, '0.1')
    textGUL.delete('1.0',tk.END)
    textGUL.insert(tk.END, '200')
    textGLL.delete('1.0',tk.END)
    textGLL.insert(tk.END, '10')
    
def run_algorithm(DB, DNA, RS, MD, SSLL, CUL, CLL, IN, MR, CR, GUL, GLL, AN, FT, ECB):
    errors = ""
    if DB == "":
        errors += 'No protein files input\n'
    if DB != "":
        file_list = DB.split("\n")
        while "" in file_list:
            file_list.remove("")
        for file in file_list:
            if isfile(file) == False:
                errors += f'{file} do not exist\n'
    if DNA == "":
        errors += 'No DNA files input\n'
    if DNA != "":
        DNA_file_list = DNA.split("\n")
        while "" in DNA_file_list:
            DNA_file_list.remove("")
        for file in DNA_file_list:
            if isfile(file) == False:
                errors += f'{file} do not exist\n'
    if RS == "":
        errors += 'No sequence patterns file inputs\n'
    if RS != "":
        RS_file_list = RS.split("\n")
        while "" in RS_file_list:
            RS_file_list.remove("")
        for file in RS_file_list:
            if isfile(file) == False:
                errors += f'{file} do not exist\n'
    if MD == "":
        errors += 'No module inputs\n'
    if SSLL == "":
        errors += 'No self score lower limit inputs\n'
    else:
        if (float(SSLL)<0) or (float(SSLL)>=1):
            errors += 'Self score lower limit should be in 0-1\n'
    if AN == "":
        errors += 'No amino acid number inputs\n'
    if FT == "":
        errors += 'No time and tempreture input\n'
    if var1.get()==0:
        if CUL == "":
            errors += 'No maximum cycles inputs\n'
        else:
            if int(CUL) < int(CLL):
                errors += 'Maximum cycles should be above minimum cycles inputs\n'
        if CLL == "":
            errors += 'No minimum cycles inputs\n'
        else:
            if int(CLL) == 0:
                errors += 'Minimum cycles should be above 0\n'
    if var1.get()==2:
        if IN == "":
            errors += 'No initial number inputs\n'
        if MR == "":
            errors += 'No mutation rate inputs\n'
        else:
            if (float(MR)<0) or (float(MR)>1):
                errors += 'Mutation rate should be in 0-1\n'
        if CR == "":
            errors += 'No crossover rate inputs\n'
        else:
            if (float(CR)<0) or (float(CR)>1):
                errors += 'Crossover rate should be in 0-1\n'
        if GUL == "":
            errors += 'No maximum generations inputs\n'
        else:
            if int(GUL) < int(GLL):
                errors += 'Maximum generations should be above minimum generations inputs\n'
        if GLL == "":
            errors += 'No minimum generations inputs\n'
        else:
            if int(GLL) == 0:
                errors += 'Minimum generations should be above 0\n'
    if errors != "":
        tk.messagebox.showerror('Error',f'{errors}')
    else:
        loadingwin = tk.Toplevel(app)
        loadingwin.title('The algorithm is running...')
        nww = 400
        nwh =80
        xaxis_nw = (sw-nww) / 2
        yaxis_nw = (sh-nwh) / 2
        loadingwin.geometry('%dx%d+%d+%d'%(nww,nwh,xaxis_nw,yaxis_nw))
        loadingbar = ttk.Progressbar(loadingwin, length=350, mode='determinate', orient=tk.HORIZONTAL)
        loadingbar.place(x=25, y=20)
        loadingbar['value'] = 0
        if var1.get()==0:
            loadingbar['maximum'] = int(CUL)+1
            Montecarlo(DB, DNA, RS, MD, int(CUL), int(CLL), float(SSLL), AN, FT, loadingwin, loadingbar, ECB)
            loadingwin.destroy()
        if var1.get()==1:
            loadingbar['maximum'] = len(MD)+1
            Greedy(DB, DNA, RS, MD, float(SSLL), AN, FT, loadingwin, loadingbar, ECB)
            loadingwin.destroy()
        if var1.get()==2:
            loadingbar['maximum'] = int(GUL)+1
            Genetic(DB, DNA, RS, MD, int(IN), float(MR), float(CR), int(GUL), int(GLL), float(SSLL), AN, FT, loadingwin, loadingbar, ECB)
            loadingwin.destroy()
        
def get_image(filename, width, height):
    im = Image.open(filename).resize((width, height))
    return ImageTk.PhotoImage(im)

def load_record():
    file = textLR.get('1.0',tk.END+"-2c")
    if isfile(file) == True:
        DBtext.delete('1.0',tk.END)
        DNAtext.delete('1.0',tk.END)
        RStext.delete('1.0',tk.END)
        records = open(file)
        FT_record=''
        tablelist = ['01h_25C','01h_37C','18h_25C','18h_37C']
        count_line = 0
        protein_line = float('inf')
        DNA_line = float('inf')
        RS_line = float('inf')
        for line in records:
            if line.startswith('This is the result from'):
                if line.split(' ')[5]=='Monte':
                    rd1.select()
                if line.split(' ')[5]=='Greedy':
                    rd2.select()
                if line.split(' ')[5]=='Genetic':
                    rd3.select()
            if line.startswith('Mostly used codons:'):
                if line.split(': ')[-1]=='ON':
                    biascheck.select()
                if line.split(': ')[-1]=='OFF':
                    biascheck.deselect()
            if line.startswith('Residue number:'):
                AN_record = line.split(': ')[-1][:-1]
                comb_jnum.current(aa_num_list.index(AN_record))
            if line.startswith('Ligation time:'):
                FT_record += line.split(': ')[-1][:-1]
                FT_record += '_'
            if line.startswith('Tempreture:'):
                FT_record += line.split(': ')[-1][:-1]
                combST.current(tablelist.index(FT_record))
            if line.startswith('Overhang quality limit:'):
                SSLL_record = line.split(': ')[-1][:-1]
                textSSLL.delete('1.0',tk.END)
                textSSLL.insert(tk.END, SSLL_record)
            if line.startswith('Maximum cycles:'):
                CUL_record = line.split(': ')[-1][:-1]
                textCUL.delete('1.0',tk.END)
                textCUL.insert(tk.END, CUL_record)
            if line.startswith('Minimum cycles:'):
                CLL_record = line.split(': ')[-1][:-1]
                textCLL.delete('1.0',tk.END)
                textCLL.insert(tk.END, CLL_record)
            if line.startswith('Maximum generations:'):
                GUL_record = line.split(': ')[-1][:-1]
                textGUL.delete('1.0',tk.END)
                textGUL.insert(tk.END, GUL_record)
            if line.startswith('Minimum generations:'):
                GLL_record = line.split(': ')[-1][:-1]
                textGLL.delete('1.0',tk.END)
                textGLL.insert(tk.END, GLL_record)
            if line.startswith('Initial number:'):
                IN_record = line.split(': ')[-1][:-1]
                textIN.delete('1.0',tk.END)
                textIN.insert(tk.END, IN_record)            
            if line.startswith('Mutation rate:'):
                MR_record = line.split(': ')[-1][:-1]
                textMR.delete('1.0',tk.END)
                textMR.insert(tk.END, MR_record)             
            if line.startswith('Crossover rate:'):
                CR_record = line.split(': ')[-1][:-1]
                textCR.delete('1.0',tk.END)
                textCR.insert(tk.END, CR_record)
            if line.startswith('Protein files:'):
                protein_line = count_line
            if (count_line>protein_line)&(line.split('.')[-1][:-1]=='fasta'):
                DBtext.insert(tk.END, line)
            if line.startswith('DNA files:'):
                protein_line = float('inf')
                DNA_line = count_line
            if (count_line>DNA_line)&(line.split('.')[-1][:-1]=='fasta'):
                DNAtext.insert(tk.END, line)
            if line.startswith('Restriction sites avoided:'):
                DNA_line = float('inf')
                RS_line = count_line
            if (count_line>RS_line)&(line.split('.')[-1][:-1]=='fasta'):
                RStext.insert(tk.END, line)
            if line.startswith('#1'):
                load_files()
                for module_record in line.split('  ')[:-1]:
                    MDbox2.insert(tk.END, module_record)    
            count_line+=1
        rightside()
        textLR.delete('1.0',tk.END)
    else:
        tk.messagebox.showerror('Error',f'{file} do not exist')

#show infos
def show_algorithm_info(event):
    win_algorithm_info.place(x=135, y=40)
def hide_algorithm_info(event):
    win_algorithm_info.place_forget()
    
def show_biascheck_info(event):
    win_biascheck_info.place(x=520, y=70)
def hide_biascheck_info(event):
    win_biascheck_info.place_forget()
    
def show_aa_info(event):
    win_aa_info.place(x=195, y=120)
def hide_aa_info(event):
    win_aa_info.place_forget()

def show_ft_info(event):
    win_ft_info.place(x=500, y=120)
def hide_ft_info(event):
    win_ft_info.place_forget()
        
def show_file_info(event):
    win_file_info.place(x=160, y=215)
def hide_file_info(event):
    win_file_info.place_forget()

def show_DNAfile_info(event):
    win_DNAfile_info.place(x=140, y=370)
def hide_DNAfile_info(event):
    win_DNAfile_info.place_forget()
    
def show_RSfile_info(event):
    win_RSfile_info.place(x=510, y=370)
def hide_RSfile_info(event):
    win_RSfile_info.place_forget()    

def show_module_info(event):
    win_module_info.place(x=130, y=490)
def hide_module_info(event):
    win_module_info.place_forget()
        
def show_AP_info(event):
    win_AP_info.place(x=450, y=40)
def hide_AP_info(event):
    win_AP_info.place_forget()       
        
def show_SSLL_info(event):
    win_SSLL_info.place(x=370, y=80)
def hide_SSLL_info(event):
    win_SSLL_info.place_forget()
        
def show_CUL_info(event):
    win_CUL_info.place(x=370, y=120)
def hide_CUL_info(event):
    win_CUL_info.place_forget()
    
def show_CLL_info(event):
    win_CLL_info.place(x=370, y=160)
def hide_CLL_info(event):
    win_CLL_info.place_forget()
        
def show_IN_info(event):
    win_IN_info.place(x=435, y=200)
def hide_IN_info(event):
    win_IN_info.place_forget()
                
def show_MR_info(event):
    win_MR_info.place(x=505, y=240)
def hide_MR_info(event):
    win_MR_info.place_forget()

def show_CR_info(event):
    win_CR_info.place(x=505, y=280)
def hide_CR_info(event):
    win_CR_info.place_forget()
                
def show_GUL_info(event):
    win_GUL_info.place(x=330, y=120)
def hide_GUL_info(event):
    win_GUL_info.place_forget()
    
def show_GLL_info(event):
    win_GLL_info.place(x=330, y=160)
def hide_GLL_info(event):
    win_GLL_info.place_forget()
    
def show_LR_info(event):
    win_LR_info.place(x=470, y=490)
def hide_LR_info(event):
    win_LR_info.place_forget()
    
# UI surface
app = tk.Tk()
app.title('Overhangs for golden gate')
sw = app.winfo_screenwidth()
sh = app.winfo_screenheight()
ww=980
wh=850
xaxis = (sw-ww) / 2-260
yaxis = (sh-wh) / 2
app.geometry('%dx%d+%d+%d'%(ww,wh,xaxis,yaxis))
app.resizable(False, False)

canvas = tk.Canvas(app, width=980, height=850)
im = get_image('appbg.jpg',980,850)
canvas.create_image(490,425,image=im)
canvas.pack()


#left side
lb1 = tk.Label(app, text='Algorithm:', fg='black', bg='PaleTurquoise', font=('Times',15,'bold'))
lb1.place(x=40,y=40,height=30,width=95)

var1 = tk.IntVar()
rd1 = tk.Radiobutton(app,text="Monte Carlo",variable=var1,value=0,font=('Times',13),borderwidth = 1, relief="raised", command=rightside)
rd1.place(x=40,y=70,height=30,width=120)
rd2 = tk.Radiobutton(app,text="Greedy",variable=var1,value=1,font=('Times',13), borderwidth = 1, relief="raised",command=rightside)
rd2.place(x=160,y=70,height=30,width=90)
rd3 = tk.Radiobutton(app,text="Genetic",variable=var1,value=2,font=('Times',13), borderwidth = 1, relief="raised",command=rightside)
rd3.place(x=250,y=70,height=30,width=90)

CheckVar = tk.IntVar()
biascheck = tk.Checkbutton(app,text='Mostly used codons',font=('Times',13),variable = CheckVar,onvalue=1,offvalue=0)
biascheck.place(x=380,y=70)
biascheck.select()

lbjnum = tk.Label(app, text='Residue number:', fg='black', bg='PaleTurquoise', font=('Times',15,'bold'))
lbjnum.place(x=40,y=120,height=30,width=155)
comb_jnum = ttk.Combobox(app,values=aa_num_list,font=('Times',12), state='readonly')
comb_jnum.place(x=40,y=150,height=30, width=200)

lbST = tk.Label(app, text='Time and tempreture:', fg='black', bg='PaleTurquoise', font=('Times',15,'bold'))
lbST.place(x=300,y=120,height=30,width=200)
combST = ttk.Combobox(app,values=['01h_25C','01h_37C','18h_25C','18h_37C'],font=('Times',12),state='readonly')
combST.place(x=300,y=150,height=30, width=200)

style1=ttk.Style() # 创建Style对象
style1.configure('1.TSeparator',background='black') # 设置背景颜色 
sep1 = ttk.Separator(app, orient='horizontal',style='1.TSeparator')  # VERTICAL为竖的分割线
sep1.place(x=40,y=200,height=1,width=350)

lb2 = tk.Label(app, text='Protein files:', fg='black', bg='MistyRose', font=('Times',15,'bold'))
lb2.place(x=40,y=215,height=30,width=120)
DBtext = tk.Text(app, bg='PowderBlue',borderwidth = 3,relief='sunken', font=('Times',13))
DBtext.place(x=40,y=245,height=100,width=400)
hook_dropfiles(DBtext , func=dragged_DBfiles)

sep2 = ttk.Separator(app, orient='horizontal',style='1.TSeparator')  # VERTICAL为竖的分割线
sep2.place(x=40,y=360,height=1,width=350)

lb4 = tk.Label(app, text='DNA files:', fg='black', bg='MistyRose', font=('Times',15,'bold'))
lb4.place(x=40,y=370,height=30,width=100)
DNAtext = tk.Text(app, bg='PowderBlue',borderwidth = 3,relief='sunken', font=('Times',13))
DNAtext.place(x=40,y=400,height=60,width=250)
hook_dropfiles(DNAtext , func=dragged_DNAfiles)

lb5 = tk.Label(app, text='Sequence patterns file:', fg='black', bg='MistyRose', font=('Times',15,'bold'))
lb5.place(x=310,y=370,height=30,width=200)
RStext = tk.Text(app, bg='PowderBlue',borderwidth = 3,relief='sunken', font=('Times',13))
RStext.place(x=310,y=400,height=60,width=250)
hook_dropfiles(RStext , func=dragged_RSfiles)

sep4 = ttk.Separator(app, orient='horizontal',style='1.TSeparator')  # VERTICAL为竖的分割线
sep4.place(x=40,y=475,height=1,width=350)

lb3 = tk.Label(app, text='Modules:', fg='black', bg='MistyRose', font=('Times',15,'bold'))
lb3.place(x=40,y=490,height=30,width=90)

SearchE = tk.Entry(app, bg='linen',borderwidth = 3, font=('Times',12))
SearchE.place(x=40, y=520,height=35, width=215)

MDbox1 = tk.Listbox(app, bg='linen',borderwidth = 3, font=('Times',12))
MDbox1.place(x=40,y=560,height=210,width=260)

MDbox2 = tk.Listbox(app, bg='linen',borderwidth = 3, font=('Times',12))
MDbox2.place(x=340,y=520,height=250,width=260)

scy1 = tk.Scrollbar(MDbox1,command=MDbox1.yview)
scy1.pack(side=tk.RIGHT, fill=tk.Y)
MDbox1.config(yscrollcommand=scy1.set)

scy2 = tk.Scrollbar(MDbox2,command=MDbox2.yview)
scy2.pack(side=tk.RIGHT, fill=tk.Y)
MDbox2.config(yscrollcommand=scy2.set)


sep3 = ttk.Separator(app, orient='vertical',style='1.TSeparator')  # VERTICAL为竖的分割线
sep3.place(x=630,y=50,height=680,width=1)

#right side
lb6 = tk.Label(app, text='Algorithm parameters:', fg='black', bg='PaleTurquoise', font=('Times',15,'bold'))
lb6.place(x=660,y=40,height=30,width=200)

lbSSLL = tk.Label(app, text='Overhang quality limit:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbSSLL.place(x=660,y=80,height=30,width=200)
textSSLL = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textSSLL.insert(tk.END, '0.1')
textSSLL.place(x=860,y=80,height=30,width=80)
#self score lower limit

lbCUL = tk.Label(app, text='Maximum cycles:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbCUL.place(x=660,y=120,height=30,width=200)
textCUL = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textCUL.insert(tk.END, '10000')
textCUL.place(x=860,y=120,height=30,width=80)
#cycle upper limit

lbCLL = tk.Label(app, text='Minimum cycles:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbCLL.place(x=660,y=160,height=30,width=200)
textCLL = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textCLL.insert(tk.END, '100')
textCLL.place(x=860,y=160,height=30,width=80)
#cycle lower limit

lbIN = tk.Label(app, text='Initial number:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbIN.place(x=660,y=200,height=30,width=200)
textIN = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textIN.insert(tk.END, '10')
textIN.place(x=860,y=200,height=30,width=80)
lbIN.place_forget()
textIN.place_forget()
#initial number

lbMR = tk.Label(app, text='Mutation rate:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbMR.place(x=660,y=240,height=30,width=200)
textMR = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textMR.insert(tk.END, '0.1')
textMR.place(x=860,y=240,height=30,width=80)
lbMR.place_forget()
textMR.place_forget()
#mutation rate

lbCR = tk.Label(app, text='Crossover rate:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbCR.place(x=660,y=280,height=30,width=200)
textCR = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textCR.insert(tk.END, '0.1')
textCR.place(x=860,y=280,height=30,width=80)
lbCR.place_forget()
textCR.place_forget()
#crossover rate

lbGUL = tk.Label(app, text='Maximum Generations:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbGUL.place(x=660,y=120,height=30,width=200)
textGUL = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textGUL.insert(tk.END, '200')
textGUL.place(x=860,y=120,height=30,width=80)
lbGUL.place_forget()
textGUL.place_forget()
#max_generation number

lbGLL = tk.Label(app, text='Minimum Generations:', fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbGLL.place(x=660,y=160,height=30,width=200)
textGLL = tk.Text(app, relief='sunken',font=('Times',15,'bold'))
textGLL.insert(tk.END, '10')
textGLL.place(x=860,y=160,height=30,width=80)
lbGLL.place_forget()
textGLL.place_forget()
#min_generation number

lbLR = tk.Label(app, text='Load report:',fg='black', bg='BlanchedAlmond', font=('Times',15,'bold'))
lbLR.place(x=660, y=490,height=30,width=130)
textLR = tk.Text(app, relief='sunken',font=('Times',13))
textLR.place(x=660, y=520, height=50, width=210)
hook_dropfiles(textLR , func=dragged_LRfiles)


#buttons
START = tk.Button(text='START', width=10, bg='SkyBlue',fg='black',font=('Times',12), command=lambda: run_algorithm(DBtext.get('1.0',tk.END+"-2c"),
                                                                                                           DNAtext.get('1.0',tk.END+"-2c"),
                                                                                                           RStext.get('1.0',tk.END+"-2c"),
                                                                                                           MDbox2.get(0,tk.END),
                                                                                                           textSSLL.get('1.0',tk.END+"-1c"),
                                                                                                           textCUL.get('1.0',tk.END+"-1c"),
                                                                                                           textCLL.get('1.0',tk.END+"-1c"),
                                                                                                           textIN.get('1.0',tk.END+"-1c"),
                                                                                                           textMR.get('1.0',tk.END+"-1c"),
                                                                                                           textCR.get('1.0',tk.END+"-1c"),
                                                                                                           textGUL.get('1.0',tk.END+"-1c"),
                                                                                                           textGLL.get('1.0',tk.END+"-1c"),
                                                                                                           comb_jnum.get(),
                                                                                                           combST.get(),
                                                                                                           CheckVar.get()
                                                                                                           ))
START.place(x=660,y=350,height=30,width=80)

RESET = tk.Button(text='RESET', bg='SkyBlue',fg='black', font=('Times',12), command=clearall)
RESET.place(x=760,y=350,height=30,width=80)

EXIT = tk.Button(text='EXIT', bg='SkyBlue',fg='black', font=('Times',12), command=app.destroy)
EXIT.place(x=860,y=350,height=30,width=80)

LOAD = tk.Button(text='LOAD', height=4, width=10, bg='SkyBlue',fg='black',font=('Times',12), command=load_files)
LOAD.place(x=450,y=250)

ADD = tk.Button(text='ADD MODULE', height=1, width=28, bg='SkyBlue',fg='black', font=('Times',12), command=add_module)
ADD.place(x=40,y=780)

DELETE = tk.Button(text='DELETE MODULE', height=1, width=16,bg='SkyBlue',fg='black', font=('Times',12), command=delete_module)
DELETE.place(x=340,y=780)

CLEAR = tk.Button(text='CLEAR ALL', height=1, width=10, bg='SkyBlue',fg='black', font=('Times',12), command=clear_module)
CLEAR.place(x=500,y=780)

SEARCH = tk.Button(text='🔍', height=1, width=3,bg='SkyBlue',fg='black', font=('Times',14), command=search_module)
SEARCH.place(x=260,y=520)

LOAD2 = tk.Button(text='LOAD', height=1, width=5, bg='SkyBlue',fg='black',font=('Times',12), command=load_record)
LOAD2.place(x=880, y=530)


#info_win
lb1.bind("<Enter>",show_algorithm_info)
lb1.bind("<Leave>",hide_algorithm_info)
win_algorithm_info = tk.Message(app,
                                text='Choose one algorithm:\nMontecarlo, Greedy, or Genetic',
                                font=('Calibri',13),
                                relief='ridge',
                                borderwidth = 3,
                                width=300,
                                bg='linen')

biascheck.bind("<Enter>",show_biascheck_info)
biascheck.bind("<Leave>",hide_biascheck_info)
win_biascheck_info = tk.Message(app,
                                text='Select this to use the most common codons as much as possible',
                                font=('Calibri',13),
                                relief='ridge',
                                borderwidth = 3,
                                width=200,
                                bg='linen')

lbjnum.bind("<Enter>",show_aa_info)
lbjnum.bind("<Leave>",hide_aa_info)
win_aa_info = tk.Message(app,
                        text='The number of residues in each junction;\neg: 1+1 means \nusing the last 1 residue in the former module\n and the first 1 residue in the latter module',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=350,
                        bg='linen')

lbST.bind("<Enter>",show_ft_info)
lbST.bind("<Leave>",hide_ft_info)
win_ft_info = tk.Message(app,
                        text='Time and tempreture selected for the ligation',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=200,
                        bg='linen')

lb2.bind("<Enter>",show_file_info)
lb2.bind("<Leave>",hide_file_info)
win_file_info = tk.Message(app,
                        text='Modules protein sequence files (fasta format)',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=200,
                        bg='linen')

lb4.bind("<Enter>",show_DNAfile_info)
lb4.bind("<Leave>",hide_DNAfile_info)
win_DNAfile_info = tk.Message(app,
                        text='Modules DNA sequence files (fasta format)',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=200,
                        bg='linen')

lb5.bind("<Enter>",show_RSfile_info)
lb5.bind("<Leave>",hide_RSfile_info)
win_RSfile_info = tk.Message(app,
                        text='Restriction sites avoided (fasta format)',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=200,
                        bg='linen')

lb3.bind("<Enter>",show_module_info)
lb3.bind("<Leave>",hide_module_info)
win_module_info = tk.Message(app,
                        text='All the modules are in the left box;\nThe modules selected are in the right box',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lb6.bind("<Enter>",show_AP_info)
lb6.bind("<Leave>",hide_AP_info)
win_AP_info = tk.Message(app,
                        text='Parameters of the selected algorithm',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=200,
                        bg='linen')

lbSSLL.bind("<Enter>",show_SSLL_info)
lbSSLL.bind("<Leave>",hide_SSLL_info)
win_SSLL_info = tk.Message(app,
                        text='discard the overhangs whose Watson−Crick pairings quality is below this percentage (0～1)',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbCUL.bind("<Enter>",show_CUL_info)
lbCUL.bind("<Leave>",hide_CUL_info)
win_CUL_info = tk.Message(app,
                        text='The maximum cycles in this algorithm;\nDefault value is 10000',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbCLL.bind("<Enter>",show_CLL_info)
lbCLL.bind("<Leave>",hide_CLL_info)
win_CLL_info = tk.Message(app,
                        text='The minimum cycles in this algorithm;\nDefault value is 100',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbIN.bind("<Enter>",show_IN_info)
lbIN.bind("<Leave>",hide_IN_info)
win_IN_info = tk.Message(app,
                        text='The initial individual number;\nDefault value is 10',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbMR.bind("<Enter>",show_MR_info)
lbMR.bind("<Leave>",hide_MR_info)
win_MR_info = tk.Message(app,
                        text='The mutation rate;\nDefault value is 0.1',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbCR.bind("<Enter>",show_CR_info)
lbCR.bind("<Leave>",hide_CR_info)
win_CR_info = tk.Message(app,
                        text='The crossover rate;\nDefault value is 0.1',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=300,
                        bg='linen')

lbGUL.bind("<Enter>",show_GUL_info)
lbGUL.bind("<Leave>",hide_GUL_info)
win_GUL_info = tk.Message(app,
                        text='The maximum generations in this algorithm;\nDefault value is 200',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=400,
                        bg='linen')

lbGLL.bind("<Enter>",show_GLL_info)
lbGLL.bind("<Leave>",hide_GLL_info)
win_GLL_info = tk.Message(app,
                        text='The minimum generations in this algorithm;\nDefault value is 10',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=400,
                        bg='linen')

lbLR.bind("<Enter>",show_LR_info)
lbLR.bind("<Leave>",hide_LR_info)
win_LR_info = tk.Message(app,
                        text='reload the arguments in a result file',
                        font=('Calibri',13),
                        relief='ridge',
                        borderwidth = 3,
                        width=180,
                        bg='linen')

tk.mainloop()
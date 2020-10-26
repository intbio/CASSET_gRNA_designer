from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re   
import random
import string
import pandas as pd
import matplotlib.pyplot as plt

import itertools

pd.set_option('display.max_columns', None)

Entrez.email = 'akgribkova@gmail.com' 

from Bio import Align
aligner = Align.PairwiseAligner()

#aligner.mode = 'global'
#aligner.open_gap_score = -1000
#aligner.match_score = 1.0
#aligner.mismatch_score = 0
#aligner.gap_score = -1000


# Definition and downloading of target genome(s)

def id_search(seq_id):
    handle=Entrez.efetch(db="nucleotide", id = seq_id, rettype='fasta')  
    for seq_record in SeqIO.parse(handle,'fasta'):
        return seq_record.seq
    
def get_pathogen_genomes(pathogen_dict):
    pathogen_genomes = {}
    for species, nc in pathogen_dict.items():
        pathogen_genomes.update( {species : id_search(nc)} )
    return(pathogen_genomes)

# Charactristics of Cas proteins

cas_info = {
    'SpCas9': {'pam_forw' : '[ACGT]GG', 'pam_rev' : 'CC[ACGT]', 'pam_len':3, 'len_target' : 20},
    'SaCas9': {'pam_forw' : '[ACGT][ACGT]G[AG][AG]T', 'pam_rev' : 'A[TC][TC]C[ACGT][ACGT]', 'pam_len':6, 'len_target' : 21}, 
    'CjCas9': {'pam_forw' : '[ACGT][ACGT][ACGT][ACGT]A[CT]AC', 'pam_rev' : 'GT[GA]T[ACGT][ACGT][ACGT][ACGT]', 'pam_len':8, 'len_target' : 22}, 
    'StCas9': { 'pam_forw' : '[ACGT][ACGT]AGAA[AT]', 'pam_rev' : '[CG]TTCT[ACGT][ACGT]', 'pam_len':7, 'len_target' : 20},
    'DpbCas12e' : { 'pam_forw' : 'TTC[ACGT]', 'pam_rev' : '[ACGT]GAA', 'pam_len':4, 'len_target' : 20},
    'AsCas12a' : { 'pam_forw' : 'TTT[ACG]', 'pam_rev' : '[TCG]AAA', 'pam_len':4, 'len_target' : 20},
}

# Distance between 2 Cas proteins in splyt-system obtained in molecular modelling

pam_out_sp_sp = pd.read_csv('spacer_lenght/PAM_out_N-N_SpCas9_complex_distances.csv')
pam_in_sp_sp = pd.read_csv('spacer_lenght/PAM_in_C-C_SpCas9_complex_distances.csv')
pam_direct_sp_sp = pd.read_csv('spacer_lenght/direct_C-N_SpCas9_complex_distances.csv')

pam_out_sp_sp = pam_out_sp_sp.loc[(pam_out_sp_sp['PAM_target independent distance']<100)&(pam_out_sp_sp['atoms']==0)&(pam_out_sp_sp['N-N_angle']<20)]['PAM_target independent distance'].values
pam_in_sp_sp = pam_in_sp_sp.loc[(pam_in_sp_sp['PAM_target independent distance']<100)&(pam_in_sp_sp['atoms']==0)&(pam_in_sp_sp['C-C_angle']<20)]['PAM_target independent distance'].values
pam_direct_sp_sp = pam_direct_sp_sp.loc[(pam_direct_sp_sp['PAM_target independent distance']<100)&(pam_direct_sp_sp['atoms']==0)&(pam_direct_sp_sp['C-N_angle']<20)]['PAM_target independent distance'].values

# Convrting distance between 2 Cas into distance between PAM starts

def system_selection_test(cas_1, cas_2, orientation, spacer = 'manual'):
    if spacer == 'avto':
        spacer_len = [19,23]
    
    if spacer == 'manual':
        if orientation =='PAM-out':
            spacer_len = pam_out_sp_sp
        elif orientation == 'PAM-in':
            spacer_len = pam_in_sp_sp
        else:
            spacer_len = pam_direct_sp_sp
            
    
    dist = []
    for i in spacer_len: #range(spacer_len[0], spacer_len[1], 1)    
        if orientation == 'PAM-out':
            dist.append( cas_info[cas_1]['pam_len'] + cas_info[cas_1]['len_target'] + i + cas_info[cas_2]['len_target'] ) 
        elif orientation == 'PAM-in' :
            dist.append( cas_info[cas_1]['pam_len'] + i ) 
        elif orientation == 'PAM_direct_forw':
            dist.append( cas_info[cas_1]['pam_len'] + i + cas_info[cas_2]['len_target'] ) 
        else: #orientation == 'PAM_direct_rev'
            dist.append( cas_info[cas_1]['pam_len'] + cas_info[cas_1]['len_target'] + i  ) 
    
    return dist

# Search PAM in forvard and reverce strand

def PAM_forw_and_rev_search(seq, pam_forw, pam_rev):
    genome = str(seq)
    PAM_positions = {}
    
    PAM_pos = []
    for m in re.finditer(pam_forw, genome):  
        PAM_pos.append(m.start())
        PAM_positions['cas_forw'] = PAM_pos
    
    PAM_pos = []  
    for m in re.finditer(pam_rev, genome): 
        PAM_pos.append(m.start())
        PAM_positions['cas_rev'] = PAM_pos           
    
    return PAM_positions

# Pathogen position for each cas and each genomes

def get_pathogen_pam(pathogen_dict, pathogen_genomes, cas_info):
    pathogen_pam = {}

    for pat in pathogen_dict.keys():
        pathogen_pam[pat] = dict(zip(cas_info.keys(), ['empty']*len(cas_info.keys())))
        for cas in cas_info.keys():
            pathogen_pam[pat][cas] =  PAM_forw_and_rev_search(pathogen_genomes[pat], cas_info[cas]['pam_forw'], cas_info[cas]['pam_rev'])
    return pathogen_pam 

# Input - pathogens, output - dict with PAM start for each pathogens 

def search_pam_pairs(pathogen_dict, pathogen_pam, cas_1, cas_2 , dist, orientation):

    if orientation == 'PAM-out':
        cas_left = 'cas_rev'
        cas_right = 'cas_forw'
    if orientation == 'PAM-in':
        cas_left = 'cas_forw'
        cas_right = 'cas_rev'
    if orientation == 'PAM-direct_forw':
        cas_left = 'cas_forw'
        cas_right = 'cas_forw'
    if orientation == 'PAM-direct_rev':
        cas_left = 'cas_rev'
        cas_right = 'cas_rev'
    
    pathogen_pam_pairs = {}
    
    for pathogen in pathogen_pam.keys():
    
        num = 0
        a = {}
        if (cas_left in pathogen_pam[pathogen][cas_1].keys()) and (cas_right in pathogen_pam[pathogen][cas_2].keys()):
            for pos_left in pathogen_pam[pathogen][cas_1][cas_left]:
                for i in dist:
                    if (pos_left + i) in pathogen_pam[pathogen][cas_2][cas_right]:
                        a[num]={pos_left: pos_left+i}
                        num+=1
            if len(a)!=0:
                pathogen_pam_pairs[pathogen] = a
    
    return pathogen_pam_pairs
        
#     
    
def search_target_pairs(pathogen_pam_pairs, pathogen_genomes, cas_1, cas_2 , orientation ) :
    
    if orientation == 'PAM-out':
        start_target_1 = cas_info[cas_1]['pam_len'] 
        end_target_1 = cas_info[cas_1]['pam_len'] + cas_info[cas_1]['len_target']
        start_target_2 = -cas_info[cas_2]['len_target'] 
        end_target_2 = 0

    if orientation == 'PAM-in':
        start_target_1 = - cas_info[cas_1]['len_target'] 
        end_target_1 = 0
        start_target_2 = cas_info[cas_2]['pam_len']
        end_target_2 = cas_info[cas_2]['pam_len'] + cas_info[cas_2]['len_target'] 
        
    if orientation == 'PAM-direct_rev':
        start_target_1 = cas_info[cas_1]['pam_len'] 
        end_target_1 = cas_info[cas_1]['pam_len'] + cas_info[cas_1]['len_target'] 
        start_target_2 = cas_info[cas_2]['pam_len'] 
        end_target_2 = cas_info[cas_2]['pam_len'] + cas_info[cas_1]['len_target']  
    
    if orientation == 'PAM-direct_forw':
        start_target_1 = - cas_info[cas_1]['len_target'] 
        end_target_1 = 0
        start_target_2 = - cas_info[cas_2]['len_target'] 
        end_target_2 = 0
    
    pathogen_targets = {}
    
    for pathogen in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)].keys():
        
        left_pos = []
        left_seq = []
        right_pos = []
        right_seq = []
        
        for d in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pathogen].keys():
            i = [*pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pathogen][d]][0]
            #for i in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pathogen][d].keys():
            left_pos.append(i)
            left_seq.append(str(pathogen_genomes[pathogen][i + start_target_1 : i + end_target_1]) )     

            k = pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pathogen][d][i]
            right_pos.append(k)
            right_seq.append(str(pathogen_genomes[pathogen][k + start_target_2 : k + end_target_2]) ) 
        pathogen_targets[pathogen] = {'left' : dict(zip(left_pos, left_seq)) }           
        pathogen_targets[pathogen]['right'] = dict(zip(right_pos, right_seq)) 
    
    return pathogen_targets

#

def get_pathogen_seq_pairs(pathogen_pam_pairs, pathogen_targets, cas_1, cas_2, orientation):
    #pathogen_seq_pairs = {}
    temp_dict={}
    for pat in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)].keys():
        #a_list = []
        #b_list = []
        pat_dict = {}
        n = 0
        for d in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pat].keys():
        #for system in pathogen_pam_pairs[pat]['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)].keys():
            for (key1, value1) in pathogen_pam_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pat][d].items():
                a = [value for (key, value) in pathogen_targets['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pat]['left'].items() if key == key1][0]
                b = [value for (key, value) in pathogen_targets['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)][pat]['right'].items() if key == value1][0]
                #a_list.append(a)
                #b_list.append(b)
                pat_dict[n] = {a:b}
                n += 1
        temp_dict[pat] = pat_dict#dict(zip(a_list, b_list))
    #pathogen_seq_pairs['{0}_{1}_{2}'.format(cas_1, cas_2, orientation)] = temp_dict
        
    return temp_dict #pathogen_seq_pairs

#

def sequence_compare(seq_a, seq_b):
    len1= len(seq_a)
    len2= len(seq_b)
    matches = 0
    for pos in range (0,min(len1,len2)) :
        if seq_a[pos] != seq_b[pos]:
            matches += 0
        else:
            matches+=1
    return matches

#

def pathogen_targets_alignment(seq1, seq2, mismatches_thres_12=0, mismatches_thres_all=0):


    alignments = aligner.align(seq1[:12], seq2[:12])
    for alignment in alignments:
        if 12 - alignment.score <= mismatches_thres_12:
            alignments = aligner.align(seq1, seq2)
            for alignment in alignments:
                if len(seq1) - alignment.score <= mismatches_thres_all:
                    return True
                
                
# на выходе индексы из pathogen_seq_pairs
def search_final_targets_left (pathogen_seq_pairs, cas_1, cas_2, orientation, list_of_pat, mismatches_thres_12=0, mismatches_thres_all=0, ): #pathogen_dict, pathogen_targets, 

    temp = {}
    i=1 #номер мишени в итоговом словаре
# итерируем по левым таргетам первого патогена
    for num in pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][list_of_pat[0]].keys():
        targets1 = [*pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][list_of_pat[0]][num] ][0]

        new_dict = {}
# переходим к следующим патогенам (но в первом все равно проверим, мб есть дупликаты таргета)
        for pat in list_of_pat:
            list_of_num = []
# итерируем по всем левым таргетам другого патогена
            for num2 in pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][pat].keys():
                targets2 = [*pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][pat][num2] ][0]
                if pathogen_targets_alignment(targets1, targets2, mismatches_thres_12=mismatches_thres_12, mismatches_thres_all=mismatches_thres_all ) == True:
                    list_of_num.append(num2)
# если для этого патогена есть совпадающие таргеты, добавляем этот лист, если нет - то не переходим к следующему патогену (сокращение времени поисков)
            if list_of_num:
                new_dict[pat] = list_of_num
            else:
                break
# если левый таргет нашелся во всех патогенах и еще не добавлен в итоговый словарь, то добавляем
        if all(pat in new_dict for pat in list_of_pat) and (new_dict not in temp.values() ):
            temp[i]=new_dict
            i+=1

    return temp 



    final_dict = {}
    k = 1

    for i in left_targets['{}_{}_{}'.format(cas_1, cas_2, orientation)].keys(): #номер мишени
        #взяли лист тартегов первого патогена 
        for target_num in left_targets['{}_{}_{}'.format(cas_1, cas_2, orientation)][i][list_of_pat[0]]: #номера одинаковых таргетов
            target1 = [*pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][list_of_pat[0]].get(target_num).values()][0]
            temp_pats = set({list_of_pat[0]})
            temp = {list_of_pat[0] : [target_num]}
            for pat in list_of_pat:
                for target_num_2 in left_targets['{}_{}_{}'.format(cas_1, cas_2, orientation)][i][pat]: #например, второй патоген, первый таргет, второй
                    target2 = [*pathogen_seq_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][pat].get(target_num_2).values()][0]
                    if pathogen_targets_alignment(target1, target2, mismatches_thres_12=mismatches_thres_12, mismatches_thres_all=mismatches_thres_all ) == True:
                        if pat not in temp.keys():
                            temp[pat] = [target_num_2]
                        else:
                            if target_num_2 not in temp[pat]:
                                temp[pat].append(target_num_2)
                        temp_pats.add(pat)
                        
            if set(list_of_pat).issubset(temp_pats):
                if temp not in final_dict.values():
                    final_dict[k] = temp
                    k+=1
    if len(final_dict)==0:
        print('no targets for {0}_{1}_{2}'.format(cas_1, cas_2, orientation))

    print(' {}_{}_{} ends: {} targets found'.format(cas_1, cas_2, orientation, len(final_dict)))
    print('_______________________')               
        
    return final_dict

def final_for_one_pat(pathogen_dict, pathogen_targets, pathogen_pam_pairs, cas_1, cas_2, orientation):
    print(' {}_{}_{} starts'.format(cas_1, cas_2, orientation))
    final = {}
    
    for target in pathogen_dict.keys():
        temp = {}
        i = 0
        print(target)
        
        if target in pathogen_pam_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)].keys():
            for k in [*pathogen_pam_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][target].keys()]:
                left_pos = [*pathogen_pam_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][target][k].keys()][0]
                right_pos = [*pathogen_pam_pairs['{}_{}_{}'.format(cas_1, cas_2, orientation)][target][k].values()][0]
                left_seq = [v for (k,v) in pathogen_targets['{}_{}_{}'.format(cas_1, cas_2, orientation)][target]['left'].items() if k ==left_pos][0]
                right_seq = [v for (k,v) in pathogen_targets['{}_{}_{}'.format(cas_1, cas_2, orientation)][target]['right'].items() if k == right_pos][0]
                temp[i] = {left_pos:left_seq, right_pos: right_seq}

                i+=1
            final[target]=temp
        
            if len(temp)==0:
                print('no targets for {0}_{1}_{2}'.format(cas_1, cas_2, orientation))

            print('in {} {}_{}_{} ends: {} targets found'.format(target, cas_1, cas_2, orientation, len(final[target])))
            print('_______________________')
        else:
            print('no targets for {0}_{1}_{2}'.format(cas_1, cas_2, orientation))
        
    return final
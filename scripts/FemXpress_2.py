# -*- coding: utf-8 -*-
#!usr/bin/python3

import sys, math
import numpy as np
import scipy.stats as stats
from sklearn.metrics.cluster import adjusted_rand_score
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
import random

inference_dir='./FemXpress/inference'
#inference_dir='/data5/wangxin/20220923fgh_mouse_brain/puplic_data/20231123_Cell_Stem_Cell_The_single-cell_and_spatial_transcriptional_landscape_of_human_gastrulation_and_early_brain_development/FemXpress/W43_subclusterX/20240220_inference/'+sys.argv[3]

dir_paths=[inference_dir]
for dir_path in dir_paths:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print('create inference dir!')
    else:
        print('inference dir already exist!')

fin = open(sys.argv[1], "r")
#f_pos_gene = open("pos_gene_map.txt", "r")
fout_1 = open("clusters_vote_method_1.tsv", "w")
fout_2 = open("clusters_vote_method_2.tsv", "w")
fout_3 = open("clusters_vote_method_3.tsv", "w")
fout_linkage = open("linkages.tsv", "w")
fout_genotype = open("genotypes.tsv", "w")
#fout_balance = open("balance.tsv", "w")
#fout_all_loci = open("all_loci.tsv", "w")
#fout_escape_1 = open("escape_gene_list_method_1.txt", "w")
#fout_escape_2 = open("escape_gene_list_method_2.txt", "w")
fout_escape_3 = open("escape_pos_list_method_3.txt", "w")

P_THRES_1 = 0.3
P_THRES_2 = 0.4
READS_THRES = 1
TOTAL_READS_AT_POSITION_THRES = 0
TOTAL_READS_COMFIRM_BASE_THRES = 0
MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND = 0
MINOR_ALLELE_PERCENTAGE_LOWER_BOUND = 0
LINKAGE_TEST_CELL_COUNT_THRES = 0
LINKAGE_TEST_PERCENTAGE_THRES = 0.8
REMOVE_SPURIOUS_LINKAGES = False
GENOTYPE_NUM_TO_RETAIN = 500
GENOTYPE_SCORE_THRES = 5
CLUSTER_DIFF = 1
PERCENTAGE_GENOTYPE_THRES = 0.02
CELL_OVERLAP_THRES = 0.5
MINIMUM_DIVIDED_CELL_COUNT = 10

PERCENTILES = []
POINT_NUM = 10
LOCAL_LINKAGE_TEST_CELL_COUNT_THRES = 0
COVERAGE_SNP_TO_STOP = [100, 30000]
LINKAGE_POSITION_NUM_THRES = 50
LINKAGE_CANDIDATE_NUM = 1000000

total_snp_count = 0
# In this dictionary, key: position on chrX, value: gene information of Gencode release 43
# Gene info: [gene_start, gene_end, gene_id, gene_type, gene_symbol]
#pos_gene_dict = {}
#for line in f_pos_gene.readlines():
#    sp = line.strip().split("\t")
#    pos_gene_dict[sp[0]] = [sp[1],sp[2],sp[3],sp[4],sp[5]]

pos_list = []
pos_stat_dict = {}
cell_dict = {}
position_dict = {}
# identification of escape positions with method one:
# if at a specific position, a cell has two bases with similar read coverage, then the cell supports the position as an escape position.
# The format of the dictionary value is [supported, not_supported]
pos_escape_1_dict = {}
# identification of escape positions with method two:
# if a position-base pair has equal chance to co-appear with another position with two different bases,
# e.g., 2A-7G/2C-7T has 50% chance to appear, 2A-7T/2C-7G also has 50% chance to appear, 
# if the case is true for a position with all other positions considered, we take the position as an escape position.
# The format of the dictionary is {pos:[num_of_supported_pairs, num_of_not_supported_pairs]}
pos_escape_2_dict = {}
# The calculated parameter lists by the percentiles of coverages of the real data
coverage_parameter_list = []
base_coverage_parameter_list = []
ratio_list = []
ari_list = []

linkage_position_num_list = []
good_linkage_dict = {}
linkage_dict = {}
position_pair_dict = {}

genotypes_list = []
scores_list = []
cell_count_dict_1 = {}
cell_count_dict_2 = {}
linkage_count_dict_1 = {}
linkage_count_dict_2 = {}
percentage_genotype_dict = {}
cell_set_1 = set()
cell_set_2 = set()

# All n (length of coverage_parameter_list) results of class 0 and class 1, each is a dictionary of cell barcodes.
# Element 0 is the reference of the results to calculate ARI. It is derived from the most rigorous coverage parameters.
method_1_all_results = []
method_2_all_results = []
method_3_all_results = []

all_escapes_method_1 = []
all_escapes_method_2 = []

all_genotypes = []
all_linkages = []


def identify_max_base(base_a, base_t, base_c, base_g):
    max_num = base_a
    max_base = "A"
    if base_t > max_num:
        max_num = base_t
        max_base = "T"
    if base_c > max_num:
        max_num = base_c
        max_base = "C"
    if base_g > max_num:
        max_num = base_g
        max_base = "G"
    return max_num, max_base

def cal_p(base_a, base_t, base_c, base_g):
    observed = [base_a, base_t, base_c, base_g]
    avg = float(base_a+base_t+base_c+base_g)/4.0
    expected = [avg, avg, avg ,avg]
    
    try:
        stat, p = stats.chisquare(f_obs=observed, f_exp=expected)
    except:
        print(observed)
        print(expected)
    return p

def cal_p_2(base_a, base_t, base_c, base_g):
    max_1=0
    max_2=0
    if base_a > max_1:
        max_1=base_a
    elif base_a > max_2:
        max_2=base_a
    if base_t > max_1:
        max_1=base_t
    elif base_t > max_2:
        max_2=base_t
    if base_c > max_1:
        max_1=base_c
    elif base_c > max_2:
        max_2=base_c
    if base_g > max_1:
        max_1=base_g
    elif base_g > max_2:
        max_2=base_g
    avg=float(max_1+max_2)/2.0
    observed = [max_1,max_2]
    expected = [avg,avg]
    stat,p = stats.chisquare(f_obs=observed, f_exp=expected)
    return p  

def find_alternative_allele(pos, base_1):
    max_cell_num = 0
    max_base = ""
    if pos in position_dict.keys():
        for base in position_dict[pos].keys():
            if base != base_1:
                if len(position_dict[pos][base]) > max_cell_num:
                    max_cell_num = len(position_dict[pos][base])
                    max_base = base
    return max_base

def find_two_alleles(pos):
    best_pos = pos
    base_1 = ""
    base_2 = ""
    if len(position_dict[pos].keys()) == 1:
        return best_pos, base_1, base_2
    else:
        max_len_base_1 = ""
        max_len_1 = 0
        max_len_base_2 = ""
        max_len_2 = 0
        for base in position_dict[pos].keys():
            if len(position_dict[pos][base]) > max_len_1:
                max_len_base_2 = max_len_base_1
                max_len_2 = max_len_1
                max_len_base_1 = base
                max_len_1 = len(position_dict[pos][base])
            elif len(position_dict[pos][base]) > max_len_2 and len(position_dict[pos][base]) <= max_len_1:
                max_len_base_2 = base
                max_len_2 = len(position_dict[pos][base])
        base_1 = max_len_base_1
        base_2 = max_len_base_2
        return best_pos, base_1, base_2

def find_best_position(candidate_position_list):
    print("Finding best starting cutting point ...")
    max_score = 0
    best_pos = candidate_position_list[0]
    base_1 = ""
    base_2 = ""
    for pos in candidate_position_list:
        if len(position_dict[pos].keys()) == 1:
            continue
        else:
            max_len_base_1 = ""
            max_len_1 = 0
            max_len_base_2 = ""
            max_len_2 = 0
            for base in position_dict[pos].keys():
                if len(position_dict[pos][base]) > max_len_1:
                    max_len_base_2 = max_len_base_1
                    max_len_2 = max_len_1
                    max_len_base_1 = base
                    max_len_1 = len(position_dict[pos][base])
                elif len(position_dict[pos][base]) > max_len_2 and len(position_dict[pos][base]) <= max_len_1:
                    max_len_base_2 = base
                    max_len_2 = len(position_dict[pos][base])
        score = (max_len_1+max_len_2) - (max_len_1-max_len_2)
        if score > max_score:
            max_score = score
            best_pos = pos
            base_1 = max_len_base_1
            base_2 = max_len_base_2
    return best_pos, base_1, base_2

def find_random_position(candidate_position_list):
    print("Finding random OK starting cutting point ...")
    base_1 = ""
    base_2 = ""
    while True:
        pos =random.choice(candidate_position_list)
        if len(position_dict[pos].keys()) == 1:
            continue
        else:
            base_1 = ""
            base_2 = ""
            for base in position_dict[pos].keys():
                if len(position_dict[pos][base]) > MINIMUM_DIVIDED_CELL_COUNT:
                    if base_1 == "":
                        base_1 = base
                    elif base_2 == "":
                        base_2 = base
                        print(pos, base_1, base_2)
                        return pos, base_1, base_2

def is_good_linkage(pos_1, base_1, pos_2, base_2):
    if pos_1 == pos_2 and base_1 == base_2:
        return True
    if pos_1 in good_linkage_dict.keys() and base_1 in good_linkage_dict[pos_1].keys() and [pos_2, base_2] in good_linkage_dict[pos_1][base_1]:
        return True
    #print("Checking ",pos_1,base_1,str(len(position_dict[pos_1][base_1])),pos_2,base_2,str(len(position_dict[pos_2][base_2])))
    alt_base_1 = find_alternative_allele(pos_1, base_1)
    alt_base_2 = find_alternative_allele(pos_2, base_2)
    #print("Alternatives ",pos_1,alt_base_1,str(len(position_dict[pos_1][alt_base_1])),pos_2,alt_base_2,str(len(position_dict[pos_2][alt_base_2])))

    pos_1_set = set()
    pos_2_set = set()
    for base in position_dict[pos_1]:
        pos_1_set.update(position_dict[pos_1][base])
        #print(pos_1, base, position_dict[pos_1][base])
    for base in position_dict[pos_2]:
        pos_2_set.update(position_dict[pos_2][base])
        #print(pos_2, base, position_dict[pos_2][base])
    intersection_1 = pos_1_set.intersection(pos_2_set)
    #print("Total cell repository: ", intersection_1)
    #for cell in intersection_1:
    #    print("Background cell: ", cell, pos_1, cell_dict[cell][pos_1], pos_2, cell_dict[cell][pos_2])

    total_cells_count = len(intersection_1)

    if total_cells_count<LOCAL_LINKAGE_TEST_CELL_COUNT_THRES:
        return False
    evident_cells_count = 0
    base_evident_cell_count=0
    alt_evident_cell_count=0
    #print("cell list: ", position_dict[pos_1][base_1])
    for cell in position_dict[pos_1][base_1]:
        #print("cell, positions, pos_2/base:", pos_2, cell_dict[cell].keys(), cell_dict[cell][pos_2])
        if pos_2 in cell_dict[cell].keys() and cell_dict[cell][pos_2] == base_2:
            evident_cells_count = evident_cells_count + 1
            base_evident_cell_count = base_evident_cell_count + 1
    #        print("Base evident in cell: ", cell, pos_1, base_1, pos_2, base_2)
    
    if alt_base_1!="" and alt_base_2!="":
        for cell in position_dict[pos_1][alt_base_1]:
            if pos_2 in cell_dict[cell].keys() and cell_dict[cell][pos_2] == alt_base_2:
                evident_cells_count = evident_cells_count + 1
                alt_evident_cell_count = alt_evident_cell_count + 1
    #            print("Alt base evident in cell: ", cell, pos_1, alt_base_1, pos_2, alt_base_2)
    #print(evident_cells_count, total_cells_count, evident_cells_count/total_cells_count)

    counter_evident_cells_count = 0
    base_counter_evident_cell_count=0
    alt_counter_evident_cell_count=0
    #print("cell list: ", position_dict[pos_1][base_1])
    if alt_base_1!="" and alt_base_2!="":
        for cell in position_dict[pos_1][alt_base_1]:
            #print("cell, positions, pos_2/base:", pos_2, cell_dict[cell].keys(), cell_dict[cell][pos_2])
            if pos_2 in cell_dict[cell].keys() and cell_dict[cell][pos_2] == base_2:
                counter_evident_cells_count = counter_evident_cells_count + 1
                base_counter_evident_cell_count = base_counter_evident_cell_count + 1
    #            print("Counter base evident in cell: ", cell, pos_1, alt_base_1, pos_2, base_2)
    
    if alt_base_1!="" and alt_base_2!="":
        for cell in position_dict[pos_1][base_1]:
            if pos_2 in cell_dict[cell].keys() and cell_dict[cell][pos_2] == alt_base_2:
                counter_evident_cells_count = counter_evident_cells_count + 1
                alt_counter_evident_cell_count = alt_counter_evident_cell_count + 1
    #            print("Counter alt base evident in cell: ", cell, pos_1, base_1, pos_2, alt_base_2)
    #print(counter_evident_cells_count, total_cells_count, counter_evident_cells_count/total_cells_count)

    #observed = [evident_cells_count/total_cells_count, counter_evident_cells_count/total_cells_count]
    #expected = [0.5, 0.5]
    #stat, p = stats.chisquare(f_obs=observed, f_exp=expected)
    #table = [[base_evident_cell_count, base_counter_evident_cell_count],[alt_evident_cell_count, alt_counter_evident_cell_count]]
    #print(table)
    #res = stats.fisher_exact(table)
    #print(res)
    #p=res[1]
    #global pos_escape_2_dict
    # if p >= P_THRES_2:
    #     if pos_1 in pos_escape_2_dict.keys():
    #         pos_escape_2_dict[pos_1][0] = pos_escape_2_dict[pos_1][0] + 1            
    #     else:
    #         pos_escape_2_dict[pos_1] = [1, 0]            
    #     if pos_2 in pos_escape_2_dict.keys():
    #         pos_escape_2_dict[pos_2][0] = pos_escape_2_dict[pos_2][0] + 1
    #     else:
    #         pos_escape_2_dict[pos_2] = [1, 0]
    # else:
    #     if pos_1 in pos_escape_2_dict.keys():
    #         pos_escape_2_dict[pos_1][1] = pos_escape_2_dict[pos_1][1] + 1
    #     else:
    #         pos_escape_2_dict[pos_1] = [0, 1] 

    #     if pos_2 in pos_escape_2_dict.keys():
    #         pos_escape_2_dict[pos_2][1] = pos_escape_2_dict[pos_2][1] + 1
    #     else:
    #         pos_escape_2_dict[pos_2] = [0, 1]

    if evident_cells_count/total_cells_count < LINKAGE_TEST_PERCENTAGE_THRES:
    #    print("Not evident!")
        return False
    
    if pos_1 in good_linkage_dict.keys():
        if base_1 in good_linkage_dict[pos_1].keys():
            good_linkage_dict[pos_1][base_1].append([pos_2, base_2])
        else:
            good_linkage_dict[pos_1][base_1] = [[pos_2,base_2]]
    else:
        good_linkage_dict[pos_1] = {}
        good_linkage_dict[pos_1][base_1] = [[pos_2,base_2]]
    #print("Good linkage: "+pos_1+","+base_1+" - "+pos_2+","+base_2)
    #fout_linkage.write(pos_1+"\t"+base_1+"\t"+pos_2+"\t"+base_2+"\n")
    return True

def is_on_the_path(pos_1, base_1, pos_2, base_2):
    if pos_2 == linkage_dict[pos_1][base_1][0][0] and base_2 == linkage_dict[pos_1][base_1][0][1]:
        return False
    pos = pos_1
    base = base_1
    while pos in linkage_dict.keys() and base in linkage_dict[pos].keys():
        if pos_2 == pos and base_2 == base:
            #print(pos_2+" "+base_2+" is already on the path.")
            return True
        next_node = linkage_dict[pos][base][0]
        pos = next_node[0]
        base = next_node[1]
    return False

def construct_genotype(best_pos, base_1, base_2, covers, selected_positions):
    #print("Constructing genotypes from "+best_pos+" "+base_1+" "+base_2+"...")
    #best_pos, base_1, base_2 = find_best_position(list(linkage_dict.keys()))
    genotypes = [[[[best_pos, base_1]], [[best_pos, base_2]]]]
    scores = [[1,1]]
    iteration_num = 0
    while True:
        updated_tag = 0
        i=0
        new_genotypes = []
        new_scores = []
        for pair in genotypes:
            if len(pair[0]) < iteration_num + 1:
                continue
            genotype_1 = pair[0]
            genotype_2 = pair[1]
            last_base_1 = genotype_1[-1]
            last_base_2 = genotype_2[-1]
            score_1 = scores[i][0]
            score_2 = scores[i][1]
            if last_base_1[0] in linkage_dict.keys() and last_base_1[1] in linkage_dict[last_base_1[0]].keys():
                for item_1 in linkage_dict[last_base_1[0]][last_base_1[1]]:
                    pos_1 = item_1[0]
                    base_1 = item_1[1]
                    if is_on_the_path(last_base_1[0], last_base_1[1], pos_1, base_1):
                        continue
                    if pos_1 in selected_positions and covers[selected_positions.index(pos_1)] == 1:
                        continue
                    if not (last_base_2[0] in linkage_dict.keys() and last_base_2[1] in linkage_dict[last_base_2[0]].keys()):
                        continue
                    for item_2 in linkage_dict[last_base_2[0]][last_base_2[1]]:
                        pos_2 = item_2[0]
                        base_2 = item_2[1]
                        if pos_1 == pos_2:
                            local_genotype_1 = genotype_1.copy()
                            local_genotype_2 = genotype_2.copy()
                            if not [pos_1, base_1] in local_genotype_1 and not [pos_2, base_2] in local_genotype_2:
                                updated_tag = 1
                                local_score_1 = score_1
                                local_score_2 = score_2
                                for local_base_1 in local_genotype_1:
                                    if [pos_1, base_1] in linkage_dict[local_base_1[0]][local_base_1[1]]:
                                        local_score_1 = local_score_1 + 1
                                for local_base_2 in local_genotype_2:
                                    if [pos_2, base_2] in linkage_dict[local_base_2[0]][local_base_2[1]]:
                                        local_score_2 = local_score_2 + 1
                                local_genotype_1.append([pos_1, base_1])
                                local_genotype_2.append([pos_2, base_2])

                                local_score_rank = 0
                                if len(new_scores) >= GENOTYPE_NUM_TO_RETAIN:
                                    for score_pair in new_scores:
                                        if score_pair[0] >= local_score_1:
                                            local_score_rank = local_score_rank + 1
                                    if local_score_rank < GENOTYPE_NUM_TO_RETAIN:
                                        if local_score_rank == 0:
                                            new_scores = [[local_score_1, local_score_2]]+new_scores[:GENOTYPE_NUM_TO_RETAIN]
                                            new_genotypes = [[local_genotype_1, local_genotype_2]]+new_genotypes[:GENOTYPE_NUM_TO_RETAIN]
                                        elif local_score_rank == GENOTYPE_NUM_TO_RETAIN-1:
                                            new_scores = new_scores[:GENOTYPE_NUM_TO_RETAIN-1]+[[local_score_1, local_score_2]]
                                            new_genotypes = new_genotypes[:GENOTYPE_NUM_TO_RETAIN-1]+[[local_genotype_1, local_genotype_2]]
                                        else:
                                            new_scores = new_scores[:local_score_rank]+[[local_score_1, local_score_2]]+new_scores[local_score_rank:GENOTYPE_NUM_TO_RETAIN]
                                        #new_scores.append([local_score_1, local_score_2])
                                            new_genotypes = new_genotypes[:local_score_rank]+[[local_genotype_1, local_genotype_2]]+new_genotypes[local_score_rank:GENOTYPE_NUM_TO_RETAIN]
                                        #new_genotypes.append([local_genotype_1, local_genotype_2])
                                else:
                                    for score_pair in new_scores:
                                        if score_pair[0] >= local_score_1:
                                            local_score_rank = local_score_rank + 1
                                    if local_score_rank == 0:
                                        new_scores = [[local_score_1, local_score_2]]+new_scores
                                        new_genotypes = [[local_genotype_1, local_genotype_2]]+new_genotypes
                                    elif local_score_rank == len(new_scores)-1:
                                        new_scores = new_scores+[[local_score_1, local_score_2]]
                                        new_genotypes = new_genotypes+[[local_genotype_1, local_genotype_2]]
                                    else:
                                        new_scores = new_scores[:local_score_rank]+[[local_score_1, local_score_2]]+new_scores[local_score_rank:]
                                        #new_scores.append([local_score_1, local_score_2])
                                        new_genotypes = new_genotypes[:local_score_rank]+[[local_genotype_1, local_genotype_2]]+new_genotypes[local_score_rank:]
                                ###print("New genotypes: "+str(local_genotype_1)+", "+str(local_genotype_2))
                                ###print("With scores: "+str(local_score_1)+", "+str(local_score_2))
                            break
            else:
                new_scores.append([score_1,score_2])
                new_genotypes.append(pair)
            i=i+1
        iteration_num = iteration_num + 1
        if updated_tag == 0:
            break
        genotypes = new_genotypes
        scores = new_scores

    ###print("Size of genotypes: "+str(len(genotypes)))
    ###print("Size of scores: "+str(len(genotypes)))
    #print("scores: "+str(scores))
    return genotypes, scores

def test_merge_two_divisions(cell_set_1, cell_set_2, local_cell_set_1, local_cell_set_2):
    intersection_1 = cell_set_1.intersection(local_cell_set_1)
    intersection_2 = cell_set_2.intersection(local_cell_set_2)
    intersection_3 = cell_set_1.intersection(local_cell_set_2)
    intersection_4 = cell_set_2.intersection(local_cell_set_1)
    c1 = (len(cell_set_1) > len(local_cell_set_1)) and len(local_cell_set_1)!=0 and (len(intersection_1)/len(local_cell_set_1)>=CELL_OVERLAP_THRES)
    c2 = (len(cell_set_1) < len(local_cell_set_1)) and len(cell_set_1)!=0 and (len(intersection_1)/len(cell_set_1)>=CELL_OVERLAP_THRES)
    c3 = (len(cell_set_2) > len(local_cell_set_2)) and len(local_cell_set_2)!=0 and (len(intersection_2)/len(local_cell_set_2)>=CELL_OVERLAP_THRES)
    c4 = (len(cell_set_2) < len(local_cell_set_2)) and len(cell_set_2)!=0 and (len(intersection_2)/len(cell_set_2)>=CELL_OVERLAP_THRES)
    c5 = (len(cell_set_1) > len(local_cell_set_2)) and len(local_cell_set_2)!=0 and (len(intersection_3)/len(local_cell_set_2)>=CELL_OVERLAP_THRES)
    c6 = (len(cell_set_1) < len(local_cell_set_2)) and len(cell_set_1)!=0 and (len(intersection_3)/len(cell_set_1)>=CELL_OVERLAP_THRES)
    c7 = (len(cell_set_2) > len(local_cell_set_1)) and len(local_cell_set_1)!=0 and (len(intersection_4)/len(local_cell_set_1)>=CELL_OVERLAP_THRES)
    c8 = (len(cell_set_2) < len(local_cell_set_1)) and len(cell_set_2)!=0 and  (len(intersection_4)/len(cell_set_2)>=CELL_OVERLAP_THRES)
    if c1 or c2 or c3 or c4:
        return 1
    elif c5 or c6 or c7 or c8:
        return -1
    else:
        return 0
    
def test_merge_list(cell_set_list):
    i = 0
    while i < len(cell_set_list)-1:
        j = i + 1
        while j < len(cell_set_list):
            if test_merge_two_divisions(cell_set_list[i][0], cell_set_list[i][1], cell_set_list[j][0], cell_set_list[j][1]) != 0:
                return True
            j = j + 1
        i = i + 1
    return False


def combine_two_divisions(cell_set_1, cell_set_2, local_cell_set_1, local_cell_set_2, cell_count_dict_1, cell_count_dict_2, local_cell_count_dict_1, local_cell_count_dict_2):
    global reverse
    intersection_1 = cell_set_1.intersection(local_cell_set_1)
    intersection_2 = cell_set_2.intersection(local_cell_set_2)
    intersection_3 = cell_set_1.intersection(local_cell_set_2)
    intersection_4 = cell_set_2.intersection(local_cell_set_1)
    c1 = (len(cell_set_1) > len(local_cell_set_1)) and len(local_cell_set_1)!=0 and (len(intersection_1)/len(local_cell_set_1)>=CELL_OVERLAP_THRES)
    c2 = (len(cell_set_1) < len(local_cell_set_1)) and len(cell_set_1)!=0 and (len(intersection_1)/len(cell_set_1)>=CELL_OVERLAP_THRES)
    c3 = (len(cell_set_2) > len(local_cell_set_2)) and len(local_cell_set_2)!=0 and (len(intersection_2)/len(local_cell_set_2)>=CELL_OVERLAP_THRES)
    c4 = (len(cell_set_2) < len(local_cell_set_2)) and len(cell_set_2)!=0 and (len(intersection_2)/len(cell_set_2)>=CELL_OVERLAP_THRES)
    c5 = (len(cell_set_1) > len(local_cell_set_2)) and len(local_cell_set_2)!=0 and (len(intersection_3)/len(local_cell_set_2)>=CELL_OVERLAP_THRES)
    c6 = (len(cell_set_1) < len(local_cell_set_2)) and len(cell_set_1)!=0 and (len(intersection_3)/len(cell_set_1)>=CELL_OVERLAP_THRES)
    c7 = (len(cell_set_2) > len(local_cell_set_1)) and len(local_cell_set_1)!=0 and (len(intersection_4)/len(local_cell_set_1)>=CELL_OVERLAP_THRES)
    c8 = (len(cell_set_2) < len(local_cell_set_1)) and len(cell_set_2)!=0 and  (len(intersection_4)/len(cell_set_2)>=CELL_OVERLAP_THRES)
    if c1 or c2 or c3 or c4:
        reverse = False
        for cell in local_cell_count_dict_1.keys():
            if cell in cell_count_dict_1.keys():
                cell_count_dict_1[cell] = cell_count_dict_1[cell] + local_cell_count_dict_1[cell]
            else:
                cell_count_dict_1[cell] = local_cell_count_dict_1[cell]
        for cell in local_cell_count_dict_2.keys():
            if cell in cell_count_dict_2.keys():
                cell_count_dict_2[cell] = cell_count_dict_2[cell] + local_cell_count_dict_2[cell]
            else:
                cell_count_dict_2[cell] = local_cell_count_dict_2[cell]
        new_cell_set_1 = set()
        new_cell_set_2 = set()
        for cell in cell_count_dict_1.keys():
            if cell in cell_count_dict_2.keys():
                if cell_count_dict_1[cell] > cell_count_dict_2[cell]:
                    new_cell_set_1.add(cell)
                elif cell_count_dict_2[cell] > cell_count_dict_1[cell]:
                    new_cell_set_2.add(cell)
            else:
                new_cell_set_1.add(cell)
        for cell in cell_count_dict_2.keys():
            if not cell in cell_count_dict_1.keys():
                new_cell_set_2.add(cell)
        return new_cell_set_1, new_cell_set_2, cell_count_dict_1, cell_count_dict_2
    elif c5 or c6 or c7 or c8:
        reverse = True
        for cell in local_cell_count_dict_1.keys():
            if cell in cell_count_dict_2.keys():
                cell_count_dict_2[cell] = cell_count_dict_2[cell] + local_cell_count_dict_1[cell]
            else:
                cell_count_dict_2[cell] = local_cell_count_dict_1[cell]
        for cell in local_cell_count_dict_2.keys():
            if cell in cell_count_dict_1.keys():
                cell_count_dict_1[cell] = cell_count_dict_1[cell] + local_cell_count_dict_2[cell]
            else:
                cell_count_dict_1[cell] = local_cell_count_dict_2[cell]
        new_cell_set_1 = set()
        new_cell_set_2 = set()
        for cell in cell_count_dict_1.keys():
            if cell in cell_count_dict_2.keys():
                if cell_count_dict_1[cell] > cell_count_dict_2[cell]:
                    new_cell_set_1.add(cell)
                elif cell_count_dict_2[cell] > cell_count_dict_1[cell]:
                    new_cell_set_2.add(cell)
            else:
                new_cell_set_1.add(cell)
        for cell in cell_count_dict_2.keys():
            if not cell in cell_count_dict_1.keys():
                new_cell_set_2.add(cell)
        return new_cell_set_1, new_cell_set_2, cell_count_dict_1, cell_count_dict_2
    return cell_set_1, cell_set_2, cell_count_dict_1, cell_count_dict_2

def merge_all_genotypes(genotype_list):
    cell_set_list = []
    cell_vote_dict_list = []
    for genotype in genotype_list:
        cell_vote_dict_1 = {}
        cell_vote_dict_2 = {}
        for pos, base in genotype[0]:
            for cell in position_dict[pos][base]:
                if cell in cell_vote_dict_1.keys():
                    cell_vote_dict_1[cell] = cell_vote_dict_1[cell] + 1
                else:
                    cell_vote_dict_1[cell] = 1
        for pos, base in genotype[1]:
            for cell in position_dict[pos][base]:
                if cell in cell_vote_dict_2.keys():
                    cell_vote_dict_2[cell] = cell_vote_dict_2[cell] + 1
                else:
                    cell_vote_dict_2[cell] = 1
        cell_set_1 = set()
        cell_set_2 = set()
        for cell in cell_vote_dict_1.keys():
            if cell in cell_vote_dict_2.keys():
                if cell_vote_dict_1[cell] > cell_vote_dict_2[cell]:
                    cell_set_1.add(cell)
                elif cell_vote_dict_2[cell] > cell_vote_dict_1[cell]:
                    cell_set_2.add(cell)
            else:
                cell_set_1.add(cell)
        for cell in cell_vote_dict_2.keys():
            if not cell in cell_vote_dict_1.keys():
                cell_set_2.add(cell)
        cell_set_list.append([cell_set_1, cell_set_2])
        cell_vote_dict_list.append([cell_vote_dict_1, cell_vote_dict_2])
    #print("cell_set_list: "+str(cell_set_list))
    print("test_merge_list: "+str(test_merge_list(cell_set_list)))
    while test_merge_list(cell_set_list):
        new_cell_set_list = cell_set_list.copy()
        new_cell_vote_dict_list = cell_vote_dict_list.copy()
        for i in range(len(new_cell_set_list)-1):
            print("Merge starts from "+str(i))
            if new_cell_set_list[i][0] != set() and new_cell_set_list[i][1] != set():
                j = i + 1
                while j < len(new_cell_set_list):
                    if new_cell_set_list[j][0] != set() and new_cell_set_list[j][1] != set():
                        if test_merge_two_divisions(new_cell_set_list[i][0], new_cell_set_list[i][1], new_cell_set_list[j][0], new_cell_set_list[j][1]) != 0:
                            print("Merge "+str(i)+" and "+str(j))
                            #print("Before merge: "+str(new_cell_set_list))
                            new_cell_set_list[i][0], new_cell_set_list[i][1], new_cell_vote_dict_list[i][0], new_cell_vote_dict_list[i][1] = combine_two_divisions(cell_set_list[i][0], cell_set_list[i][1], cell_set_list[j][0], cell_set_list[j][1], cell_vote_dict_list[i][0], cell_vote_dict_list[i][1], cell_vote_dict_list[j][0], cell_vote_dict_list[j][1])
                            new_cell_set_list[j][0] = set()
                            new_cell_set_list[j][1] = set()
                            new_cell_vote_dict_list[j][0] = {}
                            new_cell_vote_dict_list[j][1] = {}
                            #print("After merge: "+str(new_cell_set_list))
                    j = j + 1
        cell_set_list = []
        cell_vote_dict_list = []
        for i in range(len(new_cell_set_list)):
            if new_cell_set_list[i][0] != set() and new_cell_set_list[i][1] != set():
                cell_set_list.append(new_cell_set_list[i])
                cell_vote_dict_list.append(new_cell_vote_dict_list[i])
    i = 0
    largest_i = 0
    while i<len(cell_set_list):
        if len(cell_set_list[i][0]) > len(cell_set_list[largest_i][0]):
            largest_i = i
        i = i + 1
    print(largest_i)
    print(len(cell_set_list), len(cell_vote_dict_list))
    return cell_set_list[largest_i], cell_vote_dict_list[largest_i]

def find_best_divisive_position(cell_set, pos_set):
    ###print("Finding best starting divisive point ...")
    max_score = 0
    best_pos = ""
    base_1 = ""
    base_2 = ""
    for pos in pos_set:
        if len(position_dict[pos].keys()) == 1:
            continue
        else:
            max_len_base_1 = ""
            max_len_1 = 0
            max_len_base_2 = ""
            max_len_2 = 0
            for base in position_dict[pos].keys():
                all_cell_set = set(position_dict[pos][base])
                intersection = all_cell_set.intersection(cell_set)
                if len(intersection) > max_len_1:
                    max_len_base_2 = max_len_base_1
                    max_len_2 = max_len_1
                    max_len_base_1 = base
                    max_len_1 = len(intersection)
                elif len(intersection) > max_len_2 and len(intersection) <= max_len_1:
                    max_len_base_2 = base
                    max_len_2 = len(intersection)
        score = (max_len_1+max_len_2) - (max_len_1-max_len_2)
        if score > max_score:
            max_score = score
            best_pos = pos
            base_1 = max_len_base_1
            base_2 = max_len_base_2
    print(best_pos, base_1, base_2)
    return best_pos, base_1, base_2

# Define a custom function to represent your curve
def custom_curve(x, a, b, c, d):
    #return a * np.sin(b * x) + c
    return a / (1 + np.exp(-c * (x - b))) + d

# Define the derivative of the custom curve
def custom_curve_derivative(x, a, b, c, d):
    #return a * b * np.cos(b * x)
    u = c * (x - b)
    sigmoid = a / (1 + np.exp(-u)) + d
    derivative = a * c * np.exp(-u) / ((1 + np.exp(-u))**2)
    return derivative

def custom_curve_second_derivative(x, a, b, c):
    u = c * (x - b)
    return -a * c**2 * np.exp(-u) * ((1-np.exp(-u))/(1+np.exp(-u))**3)

def plotCurve(x_data, y_data):
    # try:
    #     # Fit the custom curve to the data
    #     params, covariance = curve_fit(custom_curve, x_data, y_data)
    # except:
    #     print("Cannot fit the curve to the data.")
    #     return []
    
    # Extract the fitted parameters
    # a, b, c, d = params

    # Create a smooth curve using the fitted parameters
    # x_fit = np.linspace(0, len(x_data), 1000)
    # y_fit = custom_curve(x_fit, a, b, c, d)

    # Calculate the derivative of the fitted curve
    # y_derivative = custom_curve_derivative(x_fit, a, b, c, d)
    # y_second_derivative = custom_curve_second_derivative(x_fit, a, b, c)

    # Find the approximate inflection point by locating where the derivative crosses zero
    # inflection_points = x_fit[np.where(np.diff(np.sign(y_second_derivative)))[0]]
    #print(np.sign(y_second_derivative))
    #print(np.diff(np.sign(y_second_derivative)))
    #print(np.where(np.diff(np.sign(y_second_derivative))))
    #print(x_fit[np.where(np.diff(np.sign(y_second_derivative)))[0]])

    # Plot the original data, the fitted curve, and the derivative
    plt.clf()
    plt.scatter(x_data, y_data, label="Data")
    #plt.plot(x_fit, y_fit, label="Fitted Curve", color='red')
    #plt.plot(x_fit, y_derivative, label="Derivative", color='green')
    #plt.plot(x_fit, y_second_derivative, label="Second derivative", color='blue')
    #plt.scatter(inflection_points, [custom_curve(x, a, b, c, d) for x in inflection_points], color='purple', label="Inflection Point", marker='o')
    plt.legend()
    plt.xlabel("Test rounds")
    plt.ylabel("Adjusted Rand index vs the first test")
    plt.title("Tested "+str(len(x_data))+" rounds for the optimum parameters")
    plt.savefig("inflection_v12.png")

    # print("Fitted Curve Parameters:")
    # print("a:", a)
    # print("b:", b)
    # print("c:", c)
    # print("d:", d)
    # print("Approximate Inflection Point(s):", inflection_points)

    # return inflection_points

def readData():
    global PERCENTILES, COVERAGE_SNP_TO_STOP
    global pos_list, pos_stat_dict, cell_dict, position_dict, pos_escape_1_dict, pos_escape_2_dict, total_snp_count
    pos_list = []
    pos_stat_dict = {}
    cell_dict = {}
    position_dict = {}
    # identification of escape positions with method one:
    # if at a specific position, a cell has two bases with similar read coverage, then the cell supports the position as an escape position.
    # The format of the dictionary value is [supported, not_supported]
    pos_escape_1_dict = {}
    # identification of escape positions with method two:
    # if a position-base pair has equal chance to co-appear with another position with two different bases,
    # e.g., 2A-7G/2C-7T has 50% chance to appear, 2A-7T/2C-7G also has 50% chance to appear, 
    # if the case is true for a position with all other positions considered, we take the position as an escape position.
    # The format of the dictionary is {pos:[num_of_supported_pairs, num_of_not_supported_pairs]}
    pos_escape_2_dict = {}
    i=0
    print("Preprocessing the data ...")
    for line in fin.readlines():
        if i==0:
            sp = line.strip().split("\t")
            k = 0
            for item in sp:
                #if not k==0:
                pos_list.append(item[5:])
                k=k+1
            total_snp_count = k
            print("Totally found "+str(k)+" SNPs in the sample.")
            if k<COVERAGE_SNP_TO_STOP[0]:
                break
            i=i+1
        else:
            if i%1000==0:
                print("Reached the "+str(i)+"th record")
            sp=line.strip().split("\t")
            j=0
            barcode = ""
            for item in sp:
                if j==0:
                    barcode=item
                    cell_dict[barcode] = {}
                    j=j+1
                else:
                    if item==".;.;.;.":
                        j=j+1
                        continue
                    else:
                        spsp = item.strip().split(";")
                        A_num = int(spsp[0])
                        C_num = int(spsp[1])
                        G_num = int(spsp[2])
                        T_num = int(spsp[3])
                        pos = pos_list[j-1]
                        sum = A_num + C_num + G_num + T_num
                        if pos in pos_stat_dict.keys():
                            pos_stat_dict[pos]["A"] = pos_stat_dict[pos]["A"] + A_num
                            pos_stat_dict[pos]["C"] = pos_stat_dict[pos]["C"] + C_num
                            pos_stat_dict[pos]["G"] = pos_stat_dict[pos]["G"] + G_num
                            pos_stat_dict[pos]["T"] = pos_stat_dict[pos]["T"] + T_num
                            pos_stat_dict[pos]["sum"] = pos_stat_dict[pos]["sum"] + sum
                        else:
                            pos_stat_dict[pos] = {}
                            pos_stat_dict[pos]["A"] = A_num
                            pos_stat_dict[pos]["C"] = C_num
                            pos_stat_dict[pos]["G"] = G_num
                            pos_stat_dict[pos]["T"] = T_num
                            pos_stat_dict[pos]["sum"] = sum
                        #p_val = cal_p_2(base_a=A_num, base_c=C_num, base_g=G_num, base_t=T_num)
                        #if p_val >= P_THRES_1:
                        #    if pos in pos_escape_1_dict.keys():
                        #        pos_escape_1_dict[pos][0] = pos_escape_1_dict[pos][0] + 1
                        #    else:
                        #        pos_escape_1_dict[pos] = [1, 0]
                        #else:
                        #    if pos in pos_escape_1_dict.keys():
                        #        pos_escape_1_dict[pos][1] = pos_escape_1_dict[pos][1] + 1
                        #    else:
                        #        pos_escape_1_dict[pos] = [0, 1]
                        max_num, max_base = identify_max_base(base_a=A_num, base_c=C_num, base_g=G_num, base_t=T_num)
                        #if p_val > P_THRES and max_num >= READS_THRES:
                        if True:
                            cell_dict[barcode][pos] = max_base
                            if pos in position_dict.keys():
                                if max_base in position_dict[pos].keys():
                                    position_dict[pos][max_base].append(barcode)
                                else:
                                    position_dict[pos][max_base] = [barcode]
                            else:
                                position_dict[pos] = {}
                                position_dict[pos][max_base] = [barcode]
                        j=j+1
            i=i+1
    # In all the cells of the sample, there could be at most two alleles for all cells, if there are more than two, it is wrong and we delete the minors except the two with most cells' support.
    # for pos in position_dict.keys():
    #     pos_all, base1_all, base2_all = find_two_alleles(pos)
    #     if base1_all!="" and base2_all!="":
    #         for base3 in list(position_dict[pos].keys()):
    #             if base3 != base1_all and base3 != base2_all:
    #                 for barcode in position_dict[pos][base3]:
    #                     del cell_dict[barcode][pos]
    #                 del position_dict[pos][base3]
    #         if (len(position_dict[pos][base2_all])/len(position_dict[pos][base1_all]))<MINOR_ALLELE_PERCENTAGE:
    #             for barcode in position_dict[pos][base1_all]:
    #                 del cell_dict[barcode][pos]
    #             for barcode in position_dict[pos][base2_all]:
    #                 del cell_dict[barcode][pos]
    #             del position_dict[pos][base1_all]
    #             del position_dict[pos][base2_all]

    for pos in position_dict.keys():
        pos_all, base1_all, base2_all = find_two_alleles(pos)
        if base1_all!="" and base2_all!="":
            for base3 in list(position_dict[pos].keys()):
                if base3 != base1_all and base3 != base2_all:
                    for barcode in position_dict[pos][base3]:
                        del cell_dict[barcode][pos]
                    del position_dict[pos][base3]
            ratio_list.append(len(position_dict[pos][base1_all])/len(position_dict[pos][base2_all]))

    median = np.percentile(ratio_list, 50)
    
    print("Ratio median: "+str(median))
    sns.displot(ratio_list, kde=True, bins=100)
    plt.axvline(x=median, color="red")
    plt.savefig("ratio.png")

    #for pos in position_dict.keys():
    #    pos_all, base1_all, base2_all = find_two_alleles(pos)
    #    if base1_all!="" and base2_all!="":
    #        fout_all_loci.write(pos_all+"\t"+base1_all+"\t"+str(len(position_dict[pos_all][base1_all]))+"\t"+base2_all+"\t"+str(len(position_dict[pos_all][base2_all]))+"\n")
    # print escape position list of method 1 with gene info
    #for key in sorted(list(pos_escape_1_dict.keys())):
    #    if key in pos_gene_dict.keys():
    #        fout_escape_1.write(key+"\t"+pos_gene_dict[key][4]+"\t"+pos_gene_dict[key][3]+"\t"+str(pos_escape_1_dict[key][0])+"\t"+str(pos_escape_1_dict[key][1])+"\t"+str(pos_escape_1_dict[key][0]/(pos_escape_1_dict[key][0]+pos_escape_1_dict[key][1]))+"\n")
    #    else:
    #        fout_escape_1.write(key+"\t-\t-"+"\t"+str(pos_escape_1_dict[key][0])+"\t"+str(pos_escape_1_dict[key][1])+"\t"+str(pos_escape_1_dict[key][0]/(pos_escape_1_dict[key][0]+pos_escape_1_dict[key][1]))+"\n")

    # Decide the coverage parameters according to the SNP count to retain
    coverage_list = []
    base_coverage_list = []
    for pos in pos_stat_dict.keys():
        coverage_list.append(pos_stat_dict[pos]["sum"])
        for key in pos_stat_dict[pos].keys():
            if key != "sum" and pos_stat_dict[pos][key] != 0:
                base_coverage_list.append(pos_stat_dict[pos][key])

    a = np.array(coverage_list)
    b = np.array(base_coverage_list)

    
    PERCENTILES = []
    if COVERAGE_SNP_TO_STOP[1] > total_snp_count:
        COVERAGE_SNP_TO_STOP[1] = total_snp_count
    if COVERAGE_SNP_TO_STOP[0] < total_snp_count * 0.5:
        COVERAGE_SNP_TO_STOP[0] = math.floor(total_snp_count * 0.5)
    for snp_count in range(COVERAGE_SNP_TO_STOP[0], COVERAGE_SNP_TO_STOP[1], (COVERAGE_SNP_TO_STOP[1]-COVERAGE_SNP_TO_STOP[0])//POINT_NUM):
        percentile = (1-(snp_count/total_snp_count))*100
        PERCENTILES.append(percentile)
        coverage_parameter_list.append(np.percentile(a, percentile))
        base_coverage_parameter_list.append(np.percentile(b, percentile))
    #i=0
    #for percentile in PERCENTILES:
    #    actual_count = total_snp_count*(1-(percentile/100))
    #    if actual_count>= COVERAGE_SNP_TO_STOP[0] and actual_count <= COVERAGE_SNP_TO_STOP[1]:
    #        coverage_parameter_list.append(np.percentile(a, percentile))
    #        base_coverage_parameter_list.append(np.percentile(b, percentile))
    #    i=i+1

def runPipeline(n):
    global coverage_parameter_list, base_coverage_parameter_list
    print("Total round number: "+str(len(coverage_parameter_list)))
    print("Coverage thresholds:"+str(coverage_parameter_list))
    print("Base coverage threholds: "+str(base_coverage_parameter_list))
    print("Test round "+str(n)+" ...")
    global pos_list, pos_stat_dict, cell_dict, position_dict, pos_escape_1_dict, pos_escape_2_dict
    
    global good_linkage_dict, linkage_dict, position_pair_dict
    good_linkage_dict = {}
    linkage_dict = {}
    position_pair_dict = {}

    
    global TOTAL_READS_AT_POSITION_THRES, TOTAL_READS_COMFIRM_BASE_THRES, REMOVE_SPURIOUS_LINKAGES, GENOTYPE_SCORE_THRES, MINIMUM_DIVIDED_CELL_COUNT
    if TOTAL_READS_AT_POSITION_THRES <= 0 and coverage_parameter_list != []:
        LOCAL_TOTAL_READS_AT_POSITION_THRES = coverage_parameter_list[n]
        print("TOTAL_READS_AT_POSITION_THRES not set. Automatically set to "+str(coverage_parameter_list[n])+".")

    if TOTAL_READS_COMFIRM_BASE_THRES <= 0 and base_coverage_parameter_list != []:
        LOCAL_TOTAL_READS_COMFIRM_BASE_THRES = base_coverage_parameter_list[n]
        print("TOTAL_READS_COMFIRM_BASE_THRES not set. Automatically set to "+str(base_coverage_parameter_list[n])+".")

    global ratio_list
    ratio_list = []
    for pos in position_dict.keys():
        pos_all, base1_all, base2_all = find_two_alleles(pos)
        if base1_all!="" and base2_all!="":
            valid_tag = 1
            for key in pos_stat_dict[pos].keys():
                    if key != "sum" and pos_stat_dict[pos][key] != 0:
                        if pos_stat_dict[pos][key] < LOCAL_TOTAL_READS_COMFIRM_BASE_THRES:
                            valid_tag = 0
            if pos_stat_dict[pos]["sum"] < LOCAL_TOTAL_READS_AT_POSITION_THRES:
                valid_tag = 0
            if valid_tag == 1:
                ratio_list.append(len(position_dict[pos][base1_all])/len(position_dict[pos][base2_all]))

    global MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND, MINOR_ALLELE_PERCENTAGE_LOWER_BOUND
    median = np.percentile(ratio_list, 50)
    print("Round median: "+str(median))
    diff = 45
    if MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND <= 0:
        LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND = np.percentile(ratio_list, 50+diff)
        print("MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND not set. Automatically set to "+str(LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND)+".")        
    if MINOR_ALLELE_PERCENTAGE_LOWER_BOUND <= 0:
        LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND = np.percentile(ratio_list, 50-diff)
        print("MINOR_ALLELE_PERCENTAGE_LOWER_BOUND not set. Automatically set to "+str(LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND)+".")

    global LINKAGE_TEST_CELL_COUNT_THRES, LOCAL_LINKAGE_TEST_CELL_COUNT_THRES
    if LINKAGE_TEST_CELL_COUNT_THRES <= 0:
        print("Setting background linkage parameter ...")
        total_cells_count_list = []
        for barcode in cell_dict.keys():
            positions = sorted(list(cell_dict[barcode].keys()), key=lambda x:int(x))
            i=0
            for pos_1 in positions:
                pos_check, base1_check, base2_check = find_two_alleles(pos_1)
                if base1_check!="" and base2_check!="":
                    b1_length = len(position_dict[pos_check][base1_check])
                    b2_length = len(position_dict[pos_check][base2_check])
                    b1_reads = pos_stat_dict[pos_check][base1_check]
                    b2_reads = pos_stat_dict[pos_check][base2_check]
                    condition_balanced = (b1_length >= b2_length and b1_length/b2_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b1_length/b2_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND) or (b2_length > b1_length and b2_length/b1_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b2_length/b1_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND)
                    condition_reads_confirm = b1_reads >= LOCAL_TOTAL_READS_COMFIRM_BASE_THRES and b2_reads >=LOCAL_TOTAL_READS_COMFIRM_BASE_THRES
                    if condition_balanced and condition_reads_confirm:
                        if pos_stat_dict[pos_1]["sum"] >= LOCAL_TOTAL_READS_AT_POSITION_THRES:
                            for base_1 in position_dict[pos_1].keys():
                                for pos_2 in positions[i+1:]:
                                    pos_check, base1_check, base2_check = find_two_alleles(pos_2)
                                    if base1_check!="" and base2_check!="":
                                        b1_length = len(position_dict[pos_check][base1_check])
                                        b2_length = len(position_dict[pos_check][base2_check])
                                        b1_reads = pos_stat_dict[pos_check][base1_check]
                                        b2_reads = pos_stat_dict[pos_check][base2_check]
                                        condition_balanced = (b1_length >= b2_length and b1_length/b2_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b1_length/b2_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND) or (b2_length > b1_length and b2_length/b1_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b2_length/b1_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND)
                                        condition_reads_confirm = b1_reads >= LOCAL_TOTAL_READS_COMFIRM_BASE_THRES and b2_reads >=LOCAL_TOTAL_READS_COMFIRM_BASE_THRES
                                        if condition_balanced and condition_reads_confirm:
                                            if pos_stat_dict[pos_2]["sum"] >= LOCAL_TOTAL_READS_AT_POSITION_THRES:
                                                pos_1_set = set()
                                                pos_2_set = set()
                                                for base in position_dict[pos_1]:
                                                    pos_1_set.update(position_dict[pos_1][base])
                                                for base in position_dict[pos_2]:
                                                    pos_2_set.update(position_dict[pos_2][base])
                                                intersection_1 = pos_1_set.intersection(pos_2_set)
                                                total_cells_count = len(intersection_1)
                                                total_cells_count_list.append(total_cells_count)
                i=i+1
        c = np.array(total_cells_count_list)
        print("Size of candidate linkage list: "+str(len(total_cells_count_list)))
        LOCAL_LINKAGE_CANDIDATE_NUM = LINKAGE_CANDIDATE_NUM
        if len(total_cells_count_list) < LOCAL_LINKAGE_CANDIDATE_NUM:
            LOCAL_LINKAGE_CANDIDATE_NUM = len(total_cells_count_list)
        print("LOCAL_LINKAGE_CANDIDATE_NUM set to "+str(LOCAL_LINKAGE_CANDIDATE_NUM))
        LINKAGE_PERCENTILE = (1-(LOCAL_LINKAGE_CANDIDATE_NUM/len(total_cells_count_list)))*100
        print("Linkage percentile is: "+str(LINKAGE_PERCENTILE))
        LOCAL_LINKAGE_TEST_CELL_COUNT_THRES = np.percentile(c, LINKAGE_PERCENTILE)
        if LOCAL_LINKAGE_TEST_CELL_COUNT_THRES < 4:
            LOCAL_LINKAGE_TEST_CELL_COUNT_THRES = 4
        print("LINKAGE_TEST_CELL_COUNT_THRES not set. Automatically set to: "+str(LOCAL_LINKAGE_TEST_CELL_COUNT_THRES))
    else:
        LOCAL_LINKAGE_TEST_CELL_COUNT_THRES = LINKAGE_TEST_CELL_COUNT_THRES

    print("Testing linkages ...")
    good_linkage_list = []
    j=0
    for barcode in cell_dict.keys():
        if j%1000 == 0:
            print("Tested "+str(j)+" cells for good linkages candidates ...")
        positions = sorted(list(cell_dict[barcode].keys()), key=lambda x:int(x))
        i=0
        for pos_1 in positions:
            if not pos_1 in position_pair_dict.keys():
                position_pair_dict[pos_1] = {}
            pos_check, base1_check, base2_check = find_two_alleles(pos_1)
            if base1_check!="" and base2_check!="":
                b1_length = len(position_dict[pos_check][base1_check])
                b2_length = len(position_dict[pos_check][base2_check])
                condition_balanced = (b1_length >= b2_length and b1_length/b2_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b1_length/b2_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND) or (b2_length > b1_length and b2_length/b1_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b2_length/b1_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND)
                b1_reads = pos_stat_dict[pos_check][base1_check]
                b2_reads = pos_stat_dict[pos_check][base2_check]
                condition_reads_confirm = b1_reads >= LOCAL_TOTAL_READS_COMFIRM_BASE_THRES and b2_reads >=LOCAL_TOTAL_READS_COMFIRM_BASE_THRES
                if condition_balanced and condition_reads_confirm:
                    if pos_stat_dict[pos_1]["sum"] >= LOCAL_TOTAL_READS_AT_POSITION_THRES:
                        for base_1 in position_dict[pos_1].keys():
                            if not base_1 in position_pair_dict[pos_1].keys():
                                position_pair_dict[pos_1][base_1] = {}
                            #print("Base 1: "+base_1)
                            for pos_2 in positions[i+1:]:
                                #print("Pos 2: "+pos_2)
                                if not pos_2 in position_pair_dict[pos_1][base_1].keys():
                                    position_pair_dict[pos_1][base_1][pos_2] = []
                                pos_check, base1_check, base2_check = find_two_alleles(pos_2)
                                if base1_check!="" and base2_check!="":
                                    b1_length = len(position_dict[pos_check][base1_check])
                                    b2_length = len(position_dict[pos_check][base2_check])
                                    condition_balanced = (b1_length >= b2_length and b1_length/b2_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b1_length/b2_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND) or (b2_length > b1_length and b2_length/b1_length >= LOCAL_MINOR_ALLELE_PERCENTAGE_LOWER_BOUND and b2_length/b1_length <= LOCAL_MINOR_ALLELE_PERCENTAGE_HIGHER_BOUND)
                                    b1_reads = pos_stat_dict[pos_check][base1_check]
                                    b2_reads = pos_stat_dict[pos_check][base2_check]
                                    condition_reads_confirm = b1_reads >= LOCAL_TOTAL_READS_COMFIRM_BASE_THRES and b2_reads >=LOCAL_TOTAL_READS_COMFIRM_BASE_THRES
                                    if condition_balanced and condition_reads_confirm:
                                        if pos_stat_dict[pos_2]["sum"] >= LOCAL_TOTAL_READS_AT_POSITION_THRES:
                                            for base_2 in position_dict[pos_2].keys():
                                                if base_2 in position_pair_dict[pos_1][base_1][pos_2]:
                                                    continue
                                                else:
                                                    position_pair_dict[pos_1][base_1][pos_2].append(base_2)
                                                    if is_good_linkage(pos_1, base_1, pos_2, base_2):
                                                        good_linkage_list.append([pos_1,base_1,pos_2,base_2])
                                                        barcode_cov_1 = len(position_dict[pos_1][base_1])
                                                        barcode_cov_2 = len(position_dict[pos_2][base_2])
                                                        #fout_balance.write(pos_1+"\t"+base_1+"\t"+str(barcode_cov_1)+"\t"+pos_2+"\t"+base_2+"\t"+str(barcode_cov_2)+"\n")
                                                        if pos_1 in linkage_dict.keys():
                                                            if base_1 in linkage_dict[pos_1].keys():
                                                                if not [pos_2, base_2] in linkage_dict[pos_1][base_1]:
                                                                    linkage_dict[pos_1][base_1].append([pos_2, base_2])
                                                            else:
                                                                linkage_dict[pos_1][base_1] = [[pos_2, base_2]]
                                                        else:
                                                            linkage_dict[pos_1] = {}
                                                            linkage_dict[pos_1][base_1] = [[pos_2, base_2]]
            i=i+1
        j=j+1
    all_linkages.append(good_linkage_list)

    # for key in sorted(list(pos_escape_2_dict.keys())):
    #     if key in pos_gene_dict.keys():
    #         fout_escape_2.write(key+"\t"+pos_gene_dict[key][4]+"\t"+pos_gene_dict[key][3]+"\t"+str(pos_escape_2_dict[key][0])+"\t"+str(pos_escape_2_dict[key][1])+"\t"+str(pos_escape_2_dict[key][0]/(pos_escape_2_dict[key][0]+pos_escape_2_dict[key][1]))+"\n")
    #     else:
    #         fout_escape_2.write(key+"\t-\t-"+"\t"+str(pos_escape_2_dict[key][0])+"\t"+str(pos_escape_2_dict[key][1])+"\t"+str(pos_escape_2_dict[key][0]/(pos_escape_2_dict[key][0]+pos_escape_2_dict[key][1]))+"\n")

    if REMOVE_SPURIOUS_LINKAGES:
        for pos_1 in list(linkage_dict.keys()):
            for base_1 in list(linkage_dict[pos_1].keys()):
                current_pos = linkage_dict[pos_1][base_1][0][0]
                current_index = 0
                index_list = []
                for i in range(1, len(linkage_dict[pos_1][base_1])):
                    if current_pos == linkage_dict[pos_1][base_1][i][0]:
                        if not current_index in index_list:
                            index_list.append(current_index)
                        index_list.append(i)
                    else:
                        current_pos = linkage_dict[pos_1][base_1][i][0]
                        current_index = i
                for j in sorted(index_list, reverse=True):
                    del linkage_dict[pos_1][base_1][j]
                if linkage_dict[pos_1][base_1] == []:
                    del linkage_dict[pos_1][base_1]
            if linkage_dict[pos_1] == {}:
                del linkage_dict[pos_1]

        for pos_1 in linkage_dict.keys():
            for base_1 in linkage_dict[pos_1].keys():
                for pos_2, base_2 in linkage_dict[pos_1][base_1]:
                    barcode_cov_1 = len(position_dict[pos_1][base_1])
                    barcode_cov_2 = len(position_dict[pos_2][base_2])
                    #fout_balance_after.write(pos_1+"\t"+base_1+"\t"+str(barcode_cov_1)+"\t"+pos_2+"\t"+base_2+"\t"+str(barcode_cov_2)+"\n")

    print("Get "+str(len(linkage_dict.keys()))+" positions on chrX.")
    linkage_position_num_list.append(len(linkage_dict.keys()))

    selected_positions = sorted(list(linkage_dict.keys()), key=lambda x:int(x))

    global covers, genotypes_list, scores_list, cell_count_dict_1, cell_count_dict_2, linkage_count_dict_1, linkage_count_dict_2, percentage_genotype_dict, cell_set_1, cell_set_2
    covers = [0]*len(selected_positions)
    genotypes_list = []
    scores_list = []
    cell_count_dict_1 = {}
    cell_count_dict_2 = {}
    linkage_count_dict_1 = {}
    linkage_count_dict_2 = {}
    percentage_genotype_dict = {}
    cell_set_1 = set()
    cell_set_2 = set()
    index=0
    reverse = False
    base_in_genotype = 0
    genotype_str_list=[]
    while index < len(selected_positions):
        if covers[index] == 0:
            #print("Checking "+str(index)+"th position in selected positions ...")
            pos, base_1, base_2 = find_two_alleles(selected_positions[index])
            if base_1 != "" and base_2 != "":
                genotypes, scores = construct_genotype(pos, base_1, base_2, covers, selected_positions)
                
                max_score = 0
                max_i = 0
                for i in range(len(scores)):
                    if scores[i][0] >= scores[i][1] and scores[i][0] > max_score:
                        max_score = scores[i][0]
                        max_i = i
                    elif scores[i][1] > scores[i][0] and scores[i][1] > max_score:
                        max_score = scores[i][1]
                        max_i = i

                #print("Max score: "+str(max_score))
                #print("Max i: "+str(max_i))
                #print("Max scored genotypes: "+str(genotypes[max_i]))

                for item in genotypes[max_i][0]:
                    if item[0] in selected_positions:
                        covers[selected_positions.index(item[0])] = 1
                #print("Covers: "+str(covers))
                
                if max_score > GENOTYPE_SCORE_THRES:
                    genotypes_list.append(genotypes[max_i])
                    scores_list.append(scores[max_i])
                    g_copy_0 = []
                    g_copy_1 = []
                    for locus in genotypes[max_i][0]:
                        g_copy_0.append(locus[0]+" "+locus[1])
                    for locus in genotypes[max_i][1]:
                        g_copy_1.append(locus[0]+" "+locus[1])
                    base_in_genotype = base_in_genotype + len(g_copy_0)
                    if cell_set_1 == set() and cell_set_2 == set():
                        #print("Voting for genotype 1: ")
                        for item in genotypes[max_i][0]:
                            #print(str(item)+" has voted for "+str(len(position_dict[item[0]][item[1]]))+" cells")
                            for cell in position_dict[item[0]][item[1]]:
                                if not cell in cell_count_dict_1.keys():
                                    cell_count_dict_1[cell] = 1
                                else:
                                    cell_count_dict_1[cell] = cell_count_dict_1[cell] + 1
                        #print("Voting for genotype 2: ")
                        for item in genotypes[max_i][1]:
                            #print(str(item)+" has voted for "+str(len(position_dict[item[0]][item[1]]))+" cells")
                            for cell in position_dict[item[0]][item[1]]:
                                if not cell in cell_count_dict_2.keys():
                                    cell_count_dict_2[cell] = 1
                                else:
                                    cell_count_dict_2[cell] = cell_count_dict_2[cell] + 1
                        #fout_genotype.write("|".join(g_copy_0))
                        #fout_genotype.write("&")
                        #fout_genotype.write("|".join(g_copy_1))
                        #fout_genotype.write("\n")
                        genotype_str_list.append("|".join(g_copy_0)+"&"+"|".join(g_copy_1)+"\n")
                        for barcode in cell_dict.keys():
                            index_i = 0
                            while index_i < len(g_copy_0)-1:
                                sp0 = g_copy_0[index_i].split(" ")
                                sp1 = g_copy_0[index_i+1].split(" ")
                                if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                    if barcode in linkage_count_dict_1.keys():
                                        linkage_count_dict_1[barcode] = linkage_count_dict_1[barcode] + 1
                                    else:
                                        linkage_count_dict_1[barcode] = 1
                                index_i = index_i + 1
                            index_i = 0
                            while index_i < len(g_copy_1)-1:
                                sp0 = g_copy_1[index_i].split(" ")
                                sp1 = g_copy_1[index_i+1].split(" ")
                                if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                    if barcode in linkage_count_dict_2.keys():
                                        linkage_count_dict_2[barcode] = linkage_count_dict_2[barcode] + 1
                                    else:
                                        linkage_count_dict_2[barcode] = 1
                                index_i = index_i + 1
                            for base in g_copy_0:
                                sp=base.split(" ")
                                if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                    if barcode in percentage_genotype_dict.keys():
                                        percentage_genotype_dict[barcode][0] = percentage_genotype_dict[barcode][0] + 1
                                    else:
                                        percentage_genotype_dict[barcode] = [1,0]
                            for base in g_copy_1:
                                sp=base.split(" ")
                                if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                    if barcode in percentage_genotype_dict.keys():
                                        percentage_genotype_dict[barcode][1] = percentage_genotype_dict[barcode][1] + 1
                                    else:
                                        percentage_genotype_dict[barcode] = [0,1]
                    elif cell_set_1 != set() and cell_set_2 != set():
                        #print("Combine genotypes.")
                        local_cell_count_dict_1 = {}
                        local_cell_count_dict_2 = {}
                        local_cell_set_1 = set()
                        local_cell_set_2 = set()
                        for item in genotypes[max_i][0]:
                            for cell in position_dict[item[0]][item[1]]:
                                if not cell in local_cell_count_dict_1.keys():
                                    local_cell_count_dict_1[cell] = 1
                                else:
                                    local_cell_count_dict_1[cell] = local_cell_count_dict_1[cell] + 1
                        for item in genotypes[max_i][1]:
                            for cell in position_dict[item[0]][item[1]]:
                                if not cell in local_cell_count_dict_2.keys():
                                    local_cell_count_dict_2[cell] = 1
                                else:
                                    local_cell_count_dict_2[cell] = local_cell_count_dict_2[cell] + 1
                        for cell in local_cell_count_dict_1.keys():
                            if not cell in local_cell_count_dict_2.keys():
                                local_cell_set_1.add(cell)
                            else:
                                if local_cell_count_dict_1[cell] > local_cell_count_dict_2[cell]:
                                    local_cell_set_1.add(cell)
                                elif local_cell_count_dict_1[cell] < local_cell_count_dict_2[cell]:
                                    local_cell_set_2.add(cell)
                        for cell in local_cell_count_dict_2.keys():
                            if not cell in local_cell_count_dict_1.keys():
                                local_cell_set_2.add(cell)
                        #print("Before combination: "+"genotype 1 - "+str(len(cell_count_dict_1.keys()))+" cells, genotype 2 - "+str(len(cell_count_dict_2.keys()))+" cells")
                        cell_set_1, cell_set_2, cell_count_dict_1, cell_count_dict_2 = combine_two_divisions(cell_set_1, cell_set_2, local_cell_set_1, local_cell_set_2, cell_count_dict_1, cell_count_dict_2, local_cell_count_dict_1, local_cell_count_dict_2)
                        #print("After combination: "+"genotype 1 - "+str(len(cell_count_dict_1.keys()))+" cells, genotype 2 - "+str(len(cell_count_dict_2.keys()))+" cells")
                        if reverse:
                            #fout_genotype.write("|".join(g_copy_1))
                            #fout_genotype.write("&")
                            #fout_genotype.write("|".join(g_copy_0))
                            #fout_genotype.write("\n")
                            genotype_str_list.append("|".join(g_copy_1)+"&"+"|".join(g_copy_0)+"\n")
                            for barcode in cell_dict.keys():
                                index_i = 0
                                while index_i < len(g_copy_0)-1:
                                    sp0 = g_copy_0[index_i].split(" ")
                                    sp1 = g_copy_0[index_i+1].split(" ")                            
                                    if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                        if barcode in linkage_count_dict_2.keys():
                                            linkage_count_dict_2[barcode] = linkage_count_dict_2[barcode] + 1
                                        else:
                                            linkage_count_dict_2[barcode] = 1
                                    index_i = index_i + 1
                                index_i = 0
                                while index_i < len(g_copy_1)-1:
                                    sp0 = g_copy_1[index_i].split(" ")
                                    sp1 = g_copy_1[index_i+1].split(" ")
                                    if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                        if barcode in linkage_count_dict_1.keys():
                                            linkage_count_dict_1[barcode] = linkage_count_dict_1[barcode] + 1
                                        else:
                                            linkage_count_dict_1[barcode] = 1
                                    index_i = index_i + 1
                                for base in g_copy_0:
                                    sp=base.split(" ")
                                    if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                        if barcode in percentage_genotype_dict.keys():
                                            percentage_genotype_dict[barcode][1] = percentage_genotype_dict[barcode][1] + 1
                                        else:
                                            percentage_genotype_dict[barcode] = [0,1]
                                for base in g_copy_1:
                                    sp=base.split(" ")
                                    if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                        if barcode in percentage_genotype_dict.keys():
                                            percentage_genotype_dict[barcode][0] = percentage_genotype_dict[barcode][0] + 1
                                        else:
                                            percentage_genotype_dict[barcode] = [1,0]
                        else:
                            #fout_genotype.write("|".join(g_copy_0))
                            #fout_genotype.write("&")
                            #fout_genotype.write("|".join(g_copy_1))
                            #fout_genotype.write("\n")
                            genotype_str_list.append("|".join(g_copy_0)+"&"+"|".join(g_copy_1)+"\n")
                            for barcode in cell_dict.keys():
                                index_i = 0
                                while index_i < len(g_copy_0)-1:
                                    sp0 = g_copy_0[index_i].split(" ")
                                    sp1 = g_copy_0[index_i+1].split(" ")
                                    if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                        if barcode in linkage_count_dict_1.keys():
                                            linkage_count_dict_1[barcode] = linkage_count_dict_1[barcode] + 1
                                        else:
                                            linkage_count_dict_1[barcode] = 1
                                    index_i = index_i + 1
                                index_i = 0
                                while index_i < len(g_copy_1)-1:
                                    sp0 = g_copy_1[index_i].split(" ")
                                    sp1 = g_copy_1[index_i+1].split(" ")
                                    if (sp0[0] in cell_dict[barcode].keys() and sp0[1] == cell_dict[barcode][sp0[0]]) and (sp1[0] in cell_dict[barcode].keys() and sp1[1] == cell_dict[barcode][sp1[0]]):
                                        if barcode in linkage_count_dict_2.keys():
                                            linkage_count_dict_2[barcode] = linkage_count_dict_2[barcode] + 1
                                        else:
                                            linkage_count_dict_2[barcode] = 1
                                    index_i = index_i + 1
                                for base in g_copy_0:
                                    sp=base.split(" ")
                                    if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                        if barcode in percentage_genotype_dict.keys():
                                            percentage_genotype_dict[barcode][0] = percentage_genotype_dict[barcode][0] + 1
                                        else:
                                            percentage_genotype_dict[barcode] = [1,0]
                                for base in g_copy_1:
                                    sp=base.split(" ")
                                    if sp[0] in cell_dict[barcode].keys() and cell_dict[barcode][sp[0]] == sp[1]:
                                        if barcode in percentage_genotype_dict.keys():
                                            percentage_genotype_dict[barcode][1] = percentage_genotype_dict[barcode][1] + 1
                                        else:
                                            percentage_genotype_dict[barcode] = [0,1]
                    for cell in cell_count_dict_1.keys():
                        if not cell in cell_count_dict_2.keys():
                            cell_set_1.add(cell)
                        else:
                            if cell_count_dict_1[cell] > cell_count_dict_2[cell]:
                                cell_set_1.add(cell)
                            elif cell_count_dict_1[cell] < cell_count_dict_2[cell]:
                                cell_set_2.add(cell)
                    for cell in cell_count_dict_2.keys():
                        if not cell in cell_count_dict_1.keys():
                            cell_set_2.add(cell)
        index=index+1
    if len(genotypes_list) != 0:
        cell_set, cell_count_dict = merge_all_genotypes(genotypes_list)
        cell_count_dict_1 = cell_count_dict[0]
        cell_count_dict_2 = cell_count_dict[1]
    all_genotypes.append(genotype_str_list)
    print("Size of cell count dict 1: "+str(len(cell_count_dict_1.keys())))
    print("Size of cell count dict 2: "+str(len(cell_count_dict_2.keys())))

    method_1_dict = {}
    for cell in cell_count_dict_1.keys():
        if not cell in cell_count_dict_2.keys():
            #print(cell+"\t"+str(cell_count_dict_1[cell])+"\t0\t0")
            #fout_1.write(cell+"\t"+"0"+"\n")
            if cell_count_dict_1[cell] >= CLUSTER_DIFF:
                method_1_dict[cell] = 0
        else:
            if cell_count_dict_1[cell] > cell_count_dict_2[cell] and abs(cell_count_dict_1[cell]-cell_count_dict_2[cell])>=CLUSTER_DIFF:
            #    print(cell+"\t"+str(cell_count_dict_1[cell])+"\t"+str(cell_count_dict_2[cell])+"\t0")
                #fout_1.write(cell+"\t"+"0"+"\n")
                method_1_dict[cell] = 0
            elif cell_count_dict_1[cell] < cell_count_dict_2[cell] and abs(cell_count_dict_1[cell]-cell_count_dict_2[cell])>=CLUSTER_DIFF:
            #    print(cell+"\t"+str(cell_count_dict_1[cell])+"\t"+str(cell_count_dict_2[cell])+"\t1")
                #fout_1.write(cell+"\t"+"1"+"\n")
                method_1_dict[cell] = 1
            #else:
            #    print(cell+"\t"+str(cell_count_dict_1[cell])+"\t"+str(cell_count_dict_2[cell])+"\t2")
            #    fout_1.write(cell+"\t"+"2"+"\n")

    for cell in cell_count_dict_2.keys():
        if not cell in cell_count_dict_1.keys():
            #print(cell+"\t"+"0\t"+str(cell_count_dict_2[cell])+"\t1")
            if cell_count_dict_2[cell] >= CLUSTER_DIFF:
                #fout_1.write(cell+"\t"+"1"+"\n")
                method_1_dict[cell] = 1

    method_1_all_results.append(method_1_dict)

    method_2_dict = {}
    for cell in linkage_count_dict_1.keys():
        if not cell in linkage_count_dict_2.keys():
            #print(cell+"\t"+str(linkage_count_dict_1[cell])+"\t0\t0")
            #fout_2.write(cell+"\t"+"0"+"\n")
            method_2_dict[cell] = 0
        else:
            if linkage_count_dict_1[cell] > linkage_count_dict_2[cell]:
                #print(cell+"\t"+str(linkage_count_dict_1[cell])+"\t"+str(linkage_count_dict_2[cell])+"\t0")
                #fout_2.write(cell+"\t"+"0"+"\n")
                method_2_dict[cell] = 0
            elif linkage_count_dict_1[cell] < linkage_count_dict_2[cell]:
                #print(cell+"\t"+str(linkage_count_dict_1[cell])+"\t"+str(linkage_count_dict_2[cell])+"\t1")
                #fout_2.write(cell+"\t"+"1"+"\n")
                method_2_dict[cell] = 1
            #else:
            #    print(cell+"\t"+str(linkage_count_dict_1[cell])+"\t"+str(linkage_count_dict_2[cell])+"\t2")
            #    fout_2.write(cell+"\t"+"2"+"\n")

    for cell in linkage_count_dict_2.keys():
        if not cell in linkage_count_dict_1.keys():
            #print(cell+"\t"+"0\t"+str(linkage_count_dict_2[cell])+"\t1")
            #fout_2.write(cell+"\t"+"1"+"\n")
            method_2_dict[cell] = 1

    method_2_all_results.append(method_2_dict)

    method_3_dict = {}
    for cell in percentage_genotype_dict.keys():
        if percentage_genotype_dict[cell][0] > percentage_genotype_dict[cell][1] and percentage_genotype_dict[cell][0]/base_in_genotype >= PERCENTAGE_GENOTYPE_THRES:
            #fout_3.write(cell+"\t"+"0"+"\n")
            method_3_dict[cell] = 0
        if percentage_genotype_dict[cell][0] < percentage_genotype_dict[cell][1] and percentage_genotype_dict[cell][1]/base_in_genotype >= PERCENTAGE_GENOTYPE_THRES:
            #fout_3.write(cell+"\t"+"1"+"\n")
            method_3_dict[cell] = 1
        #if percentage_genotype_dict[cell][0] == percentage_genotype_dict[cell][1]:
            #fout_3.write(cell+"\t"+"2"+"\n")
    method_3_all_results.append(method_3_dict)

    remained_cell_set = set(list(cell_dict.keys())) - cell_set_1 - cell_set_2
    remained_pos_set = set(positions) - set(selected_positions)
    best_pos, base_1, base_2 = find_best_divisive_position(remained_cell_set, remained_pos_set)
    supplement_cell_count_dict_1 = {}
    supplement_cell_count_dict_2 = {}
    supplement_cell_set_1 = set()
    supplement_cell_set_2 = set()
    while remained_cell_set != set() and remained_pos_set != set() and best_pos != "" and base_1 != "" and base_2 != "" and len(set(position_dict[best_pos][base_1]).intersection(remained_cell_set)) >= MINIMUM_DIVIDED_CELL_COUNT:
        local_supplement_cell_count_dict_1 = {}
        local_supplement_cell_count_dict_2 = {}
        local_supplement_cell_set_1 = set()
        local_supplement_cell_set_2 = set()
        local_cell_set_1 = set(position_dict[best_pos][base_1]).intersection(remained_cell_set)
        local_cell_set_2 = set(position_dict[best_pos][base_2]).intersection(remained_cell_set)
        ###print("Cell sets size: "+str(len(local_cell_set_1))+" / "+str(len(local_cell_set_2)))
        # if supplement_cell_set_1 == set() and supplement_cell_set_2 == set():
        #     for cell in local_cell_set_1:
        #         if not cell in supplement_cell_count_dict_1.keys():
        #             supplement_cell_count_dict_1[cell] = 1
        #         else:
        #             supplement_cell_count_dict_1[cell] = supplement_cell_count_dict_1[cell] + 1
        #     for cell in local_cell_set_2:
        #         if not cell in supplement_cell_count_dict_2.keys():
        #             supplement_cell_count_dict_2[cell] = 1
        #         else:
        #             supplement_cell_count_dict_2[cell] = supplement_cell_count_dict_2[cell] + 1
        # elif supplement_cell_set_1 != set() and supplement_cell_set_2 != set():
        for cell in local_cell_set_1:
            if not cell in local_supplement_cell_count_dict_1.keys():
                local_supplement_cell_count_dict_1[cell] = 1
            else:
                local_supplement_cell_count_dict_1[cell] = local_supplement_cell_count_dict_1[cell] + 1
        for cell in local_cell_set_2:
            if not cell in local_supplement_cell_count_dict_2.keys():
                local_supplement_cell_count_dict_2[cell] = 1
            else:
                local_supplement_cell_count_dict_2[cell] = local_supplement_cell_count_dict_2[cell] + 1
        for cell in local_supplement_cell_count_dict_1.keys():
            if not cell in local_supplement_cell_count_dict_2.keys():
                local_supplement_cell_set_1.add(cell)
            else:
                if local_supplement_cell_count_dict_1[cell] > local_supplement_cell_count_dict_2[cell]:
                    local_supplement_cell_set_1.add(cell)
                elif local_supplement_cell_count_dict_1[cell] < local_supplement_cell_count_dict_2[cell]:
                    local_supplement_cell_set_2.add(cell)
        for cell in local_supplement_cell_count_dict_2.keys():
            if not cell in local_supplement_cell_count_dict_1.keys():
                local_supplement_cell_set_2.add(cell)
        ###print("Intersection 1-1: ", cell_set_1.intersection(local_supplement_cell_set_1))
        ###print("Intersection 1-2: ", cell_set_1.intersection(local_supplement_cell_set_2))
        ###print("Intersection 2-1: ", cell_set_2.intersection(local_supplement_cell_set_1))
        ###print("Intersection 2-2: ", cell_set_2.intersection(local_supplement_cell_set_2))
        supplement_cell_set_1, supplement_cell_set_2, supplement_cell_count_dict_1, supplement_cell_count_dict_2 = combine_two_divisions(cell_set_1, cell_set_2, local_supplement_cell_set_1, local_supplement_cell_set_2, supplement_cell_count_dict_1, supplement_cell_count_dict_2, local_supplement_cell_count_dict_1, local_supplement_cell_count_dict_2)

        count_1 = 0
        count_2 = 0
        for cell in supplement_cell_count_dict_1.keys():
            if not cell in supplement_cell_count_dict_2.keys():
                supplement_cell_set_1.add(cell)
                count_1 = count_1 + 1
            else:
                if supplement_cell_count_dict_1[cell] > supplement_cell_count_dict_2[cell]:
                    supplement_cell_set_1.add(cell)
                    count_1 = count_1 + 1
                elif supplement_cell_count_dict_1[cell] < supplement_cell_count_dict_2[cell]:
                    supplement_cell_set_2.add(cell)
                    count_2 = count_2 + 1
        for cell in supplement_cell_count_dict_2.keys():
            if not cell in supplement_cell_count_dict_1.keys():
                supplement_cell_set_2.add(cell)
                count_2 = count_2 + 1
        #print("Currently we have "+str(count_1)+" / "+str(count_2)+" new cells.")
        remained_pos_set.remove(best_pos)
        remained_cell_set = remained_cell_set - local_cell_set_1 - local_cell_set_2
        best_pos, base_1, base_2 = find_best_divisive_position(remained_cell_set, remained_pos_set)

    for cell in supplement_cell_count_dict_1.keys():
        if not cell in supplement_cell_count_dict_2.keys():
            fout_1.write(cell+"\t"+"0"+"\n")
        else:
            if supplement_cell_count_dict_1[cell] > supplement_cell_count_dict_2[cell]:
                fout_1.write(cell+"\t"+"0"+"\n")
            elif supplement_cell_count_dict_1[cell] < supplement_cell_count_dict_2[cell]:
                fout_1.write(cell+"\t"+"1"+"\n")
            else:
                fout_1.write(cell+"\t"+"2"+"\n")

    for cell in supplement_cell_count_dict_2.keys():
        if not cell in supplement_cell_count_dict_1.keys():
            fout_1.write(cell+"\t"+"1"+"\n")

readData()
if total_snp_count<COVERAGE_SNP_TO_STOP[0]:
    print("Total SNP number is less than "+str(total_snp_count)+". The sample is not classified.")
else:
    for num_count in range(len(PERCENTILES)):
        runPipeline(num_count)
    dict_index = 0
    while dict_index < len(PERCENTILES):
        if method_1_all_results[dict_index] != {} and linkage_position_num_list[dict_index] >= LINKAGE_POSITION_NUM_THRES:
            break
        dict_index = dict_index + 1

    if dict_index == len(PERCENTILES):
        print("The sample is not classified.")
    else:
        print("Percentile "+str(PERCENTILES[dict_index])+" as reference.")
            
        j = dict_index
        while j<len(PERCENTILES):
            ref = []
            alt = []
            for cell in method_1_all_results[dict_index].keys():
                if cell in method_1_all_results[j].keys():
                    ref.append(method_1_all_results[dict_index][cell])
                    alt.append(method_1_all_results[j][cell])
            ari_list.append(adjusted_rand_score(ref, alt))
            j=j+1
    print(ari_list)
    plotCurve(range(dict_index, len(PERCENTILES)), ari_list)
    # print(inflection_point_list)
    # if len(inflection_point_list) == 0 or all(i >= 0.9 for i in ari_list):
    #     best_n = -1
    # else:
    #     best_n = math.floor(inflection_point_list[0])

    best_n = 0
    for i in range(len(ari_list)):
        if ari_list[i] >= 0.9:
            best_n = i+dict_index
        else:
            break

    print("best_n: "+str(best_n))

    for cell in method_1_all_results[best_n]:
        fout_1.write(cell+"\t"+str(method_1_all_results[best_n][cell])+"\n")
    for cell in method_2_all_results[best_n]:
        fout_2.write(cell+"\t"+str(method_2_all_results[best_n][cell])+"\n")
    for cell in method_3_all_results[best_n]:
        fout_3.write(cell+"\t"+str(method_3_all_results[best_n][cell])+"\n")

    for linkage in all_linkages[best_n]:
        fout_linkage.write(linkage[0]+"\t"+linkage[1]+"\t"+linkage[2]+"\t"+linkage[3]+"\n")

    for line in all_genotypes[best_n]:
        fout_genotype.write(line)

    
    cluster_result = method_1_all_results[best_n]
    qualified_set_1 = set()
    qualified_set_2 = set()
    for pos in pos_list:
        pos_base_count_dict_1 = {}
        pos_base_count_dict_2 = {}
        for cell in cluster_result.keys():        
            if pos in cell_dict[cell]:
                base = cell_dict[cell][pos]
                if cluster_result[cell] == 0:
                    if base in pos_base_count_dict_1.keys():
                        pos_base_count_dict_1[base] = pos_base_count_dict_1[base] + 1
                    else:
                        pos_base_count_dict_1[base] = 1
                elif cluster_result[cell] == 1:
                    if base in pos_base_count_dict_2.keys():
                        pos_base_count_dict_2[base] = pos_base_count_dict_2[base] + 1
                    else:
                        pos_base_count_dict_2[base] = 1
        if len(pos_base_count_dict_1) != 0:
            sum_count_1 = 0            
            for base in pos_base_count_dict_1.keys():
                sum_count_1 = sum_count_1 + pos_base_count_dict_1[base]
            for base in pos_base_count_dict_1.keys():
                per_1 = pos_base_count_dict_1[base]/sum_count_1
                if per_1 >= 0.2 and per_1 <= 0.8 and pos_base_count_dict_1[base] >=2:
                    qualified_set_1.add(pos)
                    break

        if len(pos_base_count_dict_2) != 0:
            sum_count_2 = 0
            for base in pos_base_count_dict_2.keys():
                sum_count_2 = sum_count_2 + pos_base_count_dict_2[base]
            for base in pos_base_count_dict_2.keys():
                per_2 = pos_base_count_dict_2[base]/sum_count_2
                if per_2 >= 0.2 and per_2 <= 0.8 and pos_base_count_dict_2[base] >=2:
                    qualified_set_2.add(pos)
                    break

    qualified_intersection = qualified_set_1.intersection(qualified_set_2)

    if len(sys.argv) == 3:
        f_annot = open(sys.argv[2], "r")
        gene_dict = {}
        for line in f_annot.readlines():
            if line[0] == "#":
                continue
            sp = line.strip().split("\t")
            if sp[0] == "chrX" and sp[2] == "gene":
                start = int(sp[3])
                end = int(sp[4])
                annot = sp[8]
                gene_dict[start] = {}
                gene_dict[start][end] = annot
        annot_pos_list = sorted(list(gene_dict.keys()))
        for pos in qualified_intersection:
            pos = int(pos)
            tag = 0
            for annot_start_pos in annot_pos_list:
                if pos >= annot_start_pos:
                    for annot_end_pos in gene_dict[annot_start_pos].keys():
                        if pos < annot_end_pos:
                            annot = gene_dict[annot_start_pos][annot_end_pos]
                            gene_name = ""
                            if "gene_name" in annot:
                                p1 = annot.index("gene_name")
                                p2 = annot[p1:].index(";")
                                
                                gene_name = annot[p1+10:p1+p2].replace("\"","")
                            elif "gene_id" in annot:
                                p1 = annot.index("gene_id")
                                p2 = annot[p1:].index(";")
                                
                                gene_name = annot[p1+8:p1+p2].replace("\"","")
                            fout_escape_3.write(str(pos)+"\t"+gene_name+"\n")
                            tag = 1
                            break
                    if tag == 1:
                        break
            if tag == 0:
                fout_escape_3.write(str(pos)+"\n")
    else:
        for pos in qualified_intersection:
            fout_escape_3.write(str(pos)+"\n")

    count_list = []
    for n in range(len(method_1_all_results)):
        fout_result = open("clusters_result_"+str(n)+".tsv", "w")
        count_0 = 0
        count_1 = 0
        # fout_test = open("test/cluster_result_round_"+str(n), "w")
        for cell in method_1_all_results[n]:
            # fout_test.write(cell+"\t"+str(method_1_all_results[n][cell])+"\n")
            if method_1_all_results[n][cell] == 0:
                count_0 = count_0 + 1
                fout_result.write(cell+"\t0\n")
            elif method_1_all_results[n][cell] == 1:
                count_1 = count_1 + 1
                fout_result.write(cell+"\t1\n")
        # fout_test.close()
        count_list.append([count_0, count_1])
        fout_result.close()
        
    
    print("Ratio of each round:")
    for item in count_list:
        if item[0] > item[1]:
            if item[1] != 0:
                print(str(item[0])+"\t"+str(item[1])+"\t"+(str(item[0]/item[1])))
            else:
                print(str(item[0])+"\t"+str(item[1])+"\t"+"0")
        else:
            if item[0] != 0:
                print(str(item[1])+"\t"+str(item[0])+"\t"+(str(item[1]/item[0])))
            else:
                print(str(item[1])+"\t"+str(item[0])+"\t"+"0")

    for round_i in range(len(all_genotypes)):
        print("Genotype of round "+str(round_i))
        for line in all_genotypes[round_i]:
            print(line.strip())

fin.close()
#f_pos_gene.close()
fout_1.close()
fout_2.close()
fout_3.close()
fout_linkage.close()
fout_genotype.close()
#fout_balance.close()
#fout_all_loci.close()
#fout_escape_1.close()
#fout_escape_2.close()
fout_escape_3.close()

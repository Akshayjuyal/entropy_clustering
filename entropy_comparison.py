'''
    File name: entropy_comparison.py
    Author: Akshay Juyal
    Date created: 05/17/2022
    Date last modified: 7/22/2022
    Python Version: 3.7.12 | packaged by conda-forge | (44db2626, Nov 15 2021, 12:53:46) [PyPy 7.3.7 with GCC 9.4.0]
'''

from collections import Counter, defaultdict
import multiprocessing as mp
from tqdm import tqdm
import random
import math
import json
from copy import deepcopy
import time
import os
import logging
from datetime import timedelta, datetime
from Bio import SeqIO
import pandas as pd
import numpy as np


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)        
    # handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def read_data(process):
    data_sequence_cluster = []
    cluster_sequence_group = defaultdict(list)

    """ pickle file with all rows having format format:
            [[seq_name, cluster_number, sequence],.....]
    """
    with open(f"/alina-data0/Akshay/entropy/data/temp_3k.pkl.npy","rb") as f:
        data_sequence_cluster = np.load(f,allow_pickle=True)

    """Grouping data with their respective cluster format:
            {cluster_number:[[seq_name, cluster_number, sequence],.....],......}
    """
    for row in data_sequence_cluster:
        cluster_sequence_group[row[1]].append([row[0],row[2]])
    cluster_sequence_group = dict(cluster_sequence_group)

    """
    cluster_info_dict: dictionary with mapping of each sequence to cluster
    size_per_cluster: dictionary with mapping of each cluster with no if sequences in it
    cluster_df_dict: python list with format: [[seq_name, cluster_number, sequence],.....]
    total_seq: total no of sequences read
    seq_len: length of each aligned sequence
    """
    # with open (f"{path}/new_clustering.json") as cl:
    #     cluster_info_dict = json.loads(cl.read())
    cluster_info_dict = dict(pd.read_pickle("/alina-data0/Akshay/entropy/data/closest_haplotypes_500k_trimmed_tn93_before_march_5_0_005"))
    size_per_cluster = Counter(cluster_info_dict.values())
    cluster_df_dict=[list(x) for x in data_sequence_cluster]
    # total_seq = len(data_sequence_cluster)
    seq_len = len(cluster_df_dict[0][2])

    return cluster_sequence_group,size_per_cluster,cluster_df_dict, seq_len

def intial_actg_per_column(seq_len,data,a_c_t_g_counts):

    for column in range (seq_len):
        # for seq_name in fasta_df["sequence_name"]:
        for row in data:
            cluster= row[1]
            dict_len = len(a_c_t_g_counts[cluster])
            if row[2][column]== "-":
                if  column < dict_len:
                    a_c_t_g_counts[cluster][column]["-"]+=1
                else:
                    a_c_t_g_counts[cluster].append({"A":0,"C":0,"T":0,"G":0,"-":1,"h_dist":0,"entropy":0})
            elif row[2][column]== "A":
                if column < dict_len:
                    a_c_t_g_counts[cluster][column]["A"]+=1
                else:
                    a_c_t_g_counts[cluster].append({"A":1,"C":0,"T":0,"G":0,"-":0,"h_dist":0,"entropy":0})
            elif row[2][column]== "C":
                if column < dict_len:
                    a_c_t_g_counts[cluster][column]["C"]+=1
                else:
                    a_c_t_g_counts[cluster].append({"A":0,"C":1,"T":0,"G":0,"-":0,"h_dist":0,"entropy":0})
            elif row[2][column]== "T":
                if column < dict_len:
                    a_c_t_g_counts[cluster][column]["T"]+=1
                else:
                    a_c_t_g_counts[cluster].append({"A":0,"C":0,"T":1,"G":0,"-":0,"h_dist":0,"entropy":0})
            elif row[2][column]== "G":
                if column < dict_len:
                    a_c_t_g_counts[cluster][column]["G"]+=1
                else:
                    a_c_t_g_counts[cluster].append({"A":0,"C":0,"T":0,"G":1,"-":0,"h_dist":0,"entropy":0})
    return a_c_t_g_counts
            
def calc_initial_hamming_dist_and_entropy_per_column(seq_len,size_per_cluster,a_c_t_g_counts):
    # import pdb;pdb.set_trace()
    
    for column in range (seq_len):
        
        for cluster in [*size_per_cluster]:
            # print((cluster))
            entropy = 0
            max_actg = max([x[1] for x in a_c_t_g_counts[cluster][column].items() if x[0] in ["A","C","T","G"]])
            dash = a_c_t_g_counts[cluster][column]["-"]
            a_c_t_g_counts[cluster][column]["h_dist"] = size_per_cluster[cluster]-(max_actg+dash)
            # hamming = hamming + a_c_t_g_counts[cluster][column]["h_dist"]
            
            total = sum([x[1] for x in a_c_t_g_counts[cluster][column].items() if x[0] in ["A","C","T","G"]])
            if (total > 0):
                n_a_normalized = a_c_t_g_counts[cluster][column]["A"] / total
                #print("n_a_normalized =" , n_a_normalized)
                n_c_normalized = a_c_t_g_counts[cluster][column]["C"] / total
                n_t_normalized = a_c_t_g_counts[cluster][column]["T"] / total
                n_g_normalized = a_c_t_g_counts[cluster][column]["G"] / total
                entropy_a = n_a_normalized * math.log2(n_a_normalized) if n_a_normalized > 0 else 0
                entropy_c = n_c_normalized * math.log2(n_c_normalized) if n_c_normalized > 0 else 0
                entropy_g = n_t_normalized * math.log2(n_t_normalized) if n_t_normalized > 0 else 0
                entropy_t = n_g_normalized * math.log2(n_g_normalized) if n_g_normalized > 0 else 0
                entropy = entropy + entropy_a + entropy_c + entropy_g + entropy_t
            try:
                a_c_t_g_counts[cluster][column]["entropy"] = abs(entropy)
            except:
                print(entropy)
                break
        # hamming_entropy_dict[cluster]={
        #     "hamming":hamming,
        #     "entropy":entropy
        # }
    return a_c_t_g_counts

def calc_initial_hamming_entropy(size_per_cluster, seq_len, a_c_t_g_counts):
    hamming_entropy_dict_per_cluster={}
    for cluster_id in [*size_per_cluster]:
        entropy, hamming, total_seq= 0, 0,sum(size_per_cluster[key] for key in size_per_cluster.keys())
        for val in range(seq_len):
            hamming = hamming + a_c_t_g_counts[cluster_id][val]["h_dist"]
            entropy = entropy + a_c_t_g_counts[cluster_id][val]["entropy"]
            if hamming < 0:
                print(cluster_id,val) 
                break
        hamming_entropy_dict_per_cluster[cluster_id] = {
            "hamming":hamming*((size_per_cluster[cluster_id]/total_seq)),
            "entropy":entropy
        }
    return hamming_entropy_dict_per_cluster

def get_initial_calculations(process):
    cluster_sequence_group,size_per_cluster,cluster_df_dict, seq_len = read_data(process)

    init_a_c_t_g_counts =  dict(zip(list(size_per_cluster.keys()),[[] for val in range(len(size_per_cluster.keys()))]))

    intrim_a_c_t_g_counts = intial_actg_per_column(seq_len,cluster_df_dict,init_a_c_t_g_counts)
    a_c_t_g_counts = calc_initial_hamming_dist_and_entropy_per_column(seq_len,size_per_cluster,intrim_a_c_t_g_counts)
    hamming_entropy_dict_per_cluster = calc_initial_hamming_entropy(size_per_cluster, seq_len, a_c_t_g_counts)

    return hamming_entropy_dict_per_cluster,cluster_sequence_group,a_c_t_g_counts, size_per_cluster

# swapping sequences 

def compare_entropy_on_swap(old_size_per_cluster, size_per_cluster, new_entropy_per_cluster,old_entropy_per_cluster):
    new_en, old_en,total_seq = 0, 0 ,sum(size_per_cluster[key] for key in size_per_cluster.keys())
    for cluster_id in [*size_per_cluster]:
        new_en +=new_entropy_per_cluster[cluster_id]["entropy"]*((size_per_cluster[cluster_id]/total_seq))
        old_en +=old_entropy_per_cluster[cluster_id]["entropy"] *(( old_size_per_cluster[cluster_id]/total_seq))

    return new_en, old_en

def cal_swap_hamming_entropy_per_cluster(seq_len,temp_a_c_t_g_counts,size_per_cluster):
    temp_calculations_hamming_entropy={}
    for cluster_id in [*size_per_cluster]:
        entropy, hamming, total_seq= 0, 0,sum(size_per_cluster[key] for key in size_per_cluster.keys())
        for val in range(seq_len):
            hamming = hamming + temp_a_c_t_g_counts[cluster_id][val]["h_dist"]
            entropy = entropy + temp_a_c_t_g_counts[cluster_id][val]["entropy"]
            if hamming < 0:
                print("======>",cluster_id,val)
        
        temp_calculations_hamming_entropy[cluster_id] = {
            "hamming":hamming*((size_per_cluster[cluster_id]/total_seq)) if size_per_cluster[cluster_id] != 0 else 0,
            "entropy":entropy
        }
    return temp_calculations_hamming_entropy

def calc_swap_hamming_dist_and_entropy_per_column(seq_len,temp_a_c_t_g_counts, first_r_cluster, second_r_cluster,size_per_cluster):
    
    for column in range (seq_len):
        
        for cluster in [first_r_cluster, second_r_cluster]:
            # print((cluster))
            entropy = 0
            max_actg = max([x[1] for x in temp_a_c_t_g_counts[cluster][column].items() if x[0] in ["A","C","T","G"]])
            dash = temp_a_c_t_g_counts[cluster][column]["-"]
            temp_a_c_t_g_counts[cluster][column]["h_dist"] = size_per_cluster[cluster]-(max_actg+dash)
            # hamming = hamming + temp_a_c_t_g_counts[cluster][column]["h_dist"]
            
            total = sum([x[1] for x in temp_a_c_t_g_counts[cluster][column].items() if x[0] in ["A","C","T","G"]])
            if (total > 0):
                n_a_normalized = temp_a_c_t_g_counts[cluster][column]["A"] / total
                #print("n_a_normalized =" , n_a_normalized)
                n_c_normalized = temp_a_c_t_g_counts[cluster][column]["C"] / total
                n_t_normalized = temp_a_c_t_g_counts[cluster][column]["T"] / total
                n_g_normalized = temp_a_c_t_g_counts[cluster][column]["G"] / total
                entropy_a = n_a_normalized * math.log2(n_a_normalized) if n_a_normalized > 0 else 0
                entropy_c = n_c_normalized * math.log2(n_c_normalized) if n_c_normalized > 0 else 0
                entropy_g = n_t_normalized * math.log2(n_t_normalized) if n_t_normalized > 0 else 0
                entropy_t = n_g_normalized * math.log2(n_g_normalized) if n_g_normalized > 0 else 0
                entropy = entropy + entropy_a + entropy_c + entropy_g + entropy_t
            temp_a_c_t_g_counts[cluster][column]["entropy"] = abs(entropy)
            
    return temp_a_c_t_g_counts

def recount_actg_after_removing(temp_a_c_t_g_counts,sequence, first):
    for idx, val in enumerate(sequence):
        if val== "-":
            temp_a_c_t_g_counts[first][idx]["-"]-=1
        elif val== "A":
            temp_a_c_t_g_counts[first][idx]["A"]-=1
        elif val== "C":
            temp_a_c_t_g_counts[first][idx]["C"]-=1
        elif val== "T":
            temp_a_c_t_g_counts[first][idx]["T"]-=1
        elif val== "G":
            temp_a_c_t_g_counts[first][idx]["G"]-=1
    
    return temp_a_c_t_g_counts

def recount_actg_after_adding(temp_a_c_t_g_counts,sequence, second):
    for idx, val in enumerate(sequence):
        if val== "-":
            temp_a_c_t_g_counts[second][idx]["-"]+=1
        elif val== "A":
            temp_a_c_t_g_counts[second][idx]["A"]+=1
        elif val== "C":
            temp_a_c_t_g_counts[second][idx]["C"]+=1
        elif val== "T":
            temp_a_c_t_g_counts[second][idx]["T"]+=1
        elif val== "G":
            temp_a_c_t_g_counts[second][idx]["G"]+=1
    
    return temp_a_c_t_g_counts

def pick_random_seq_and_cluster(cluster_sequence_group):
    """Picking Random clusters to swap sequences"""
    all_clusters = list(key for key in cluster_sequence_group.keys() if len(cluster_sequence_group[key]) )
    first_r_cluster = random.choice(all_clusters)
    all_clusters.remove(first_r_cluster)
    second_r_cluster = random.choice(all_clusters)
    """Picking a random sequnce from first random cluster"""
    random_seq_to_swap = random.choice(range(len(cluster_sequence_group[first_r_cluster])))

    return first_r_cluster,second_r_cluster,random_seq_to_swap

def final_funtion(iteration_info_dict,
                    a_c_t_g_counts,orignal_a_c_t_g_counts,
                    orignal_clusters, cluster_sequence_group,
                    orignal_calculations_hamming_entropy, hamming_entropy_dict_per_cluster, process, snapshot_counter = -20):
    
    orignal_clust = {y[0]:x for x in [*orignal_clusters] for y in orignal_clusters[x]}
    new_clust = {y[0]:x for x in [*cluster_sequence_group] for y in cluster_sequence_group[x]}
    if not snapshot_counter ==-20:
        path = f"/alina-data0/Akshay/entropy/data/monte_carlo/tags/1k_pert/run1/{snapshot_counter}"
        if not os.path.exists(path):
            os.mkdir(path)  
    else:
        path = f"/alina-data0/Akshay/entropy/data/comparison_tags_with_time/24hour/30k/{process}"
        if not os.path.exists(path):
            os.mkdir(path)


    with open (f"{path}/iteration_info.json","w") as iid:
        json.dump(iteration_info_dict,iid)
    pd.DataFrame.from_dict(iteration_info_dict, orient='index').to_csv(f"{path}/iteration_info.csv",index=False)
    
    with open (f"{path}/new_actg_counts.json","w") as actg:
        json.dump(a_c_t_g_counts,actg)
    with open (f"{path}/new_clustering.json","w") as nc:
        json.dump(new_clust,nc)
    
    with open (f"{path}/new_ent_ham_per_cluster.json","w") as neh:
        json.dump(hamming_entropy_dict_per_cluster,neh)
    
    if snapshot_counter ==-20:
        with open (f"{path}/orig_actg_counts.json","w") as oactg:
            json.dump(orignal_a_c_t_g_counts,oactg)

        with open (f"{path}/orignal_clustering.json","w") as oc:
            json.dump(orignal_clust,oc)

        with open (f"{path}/orignal_ent_ham_per_cluster.json","w") as oeh:
            json.dump(orignal_calculations_hamming_entropy,oeh)


def calculate_entropy_on_swap(process):
    """
    calculated once and stored

    orignal_calculations_hamming_entropy 
    orignal_clusters
    orignal_a_c_t_g_counts

    used for comparing swap entropy

    temp_calculations_hamming_entropy
    temp_cluster_sequence
    temp_a_c_t_g_counts

    used for saving final data if entropy improves

    hamming_entropy_dict_per_cluster
    cluster_sequence_group
    a_c_t_g_counts

    """
    logfile = f"/alina-data0/Akshay/entropy/data/comparison_tags_with_time/24hour/30k/{process}/logname_{process}.log"
    logger = setup_logger(f"logname_{process}", logfile, logging.DEBUG)

    hamming_entropy_dict_per_cluster,cluster_sequence_group,\
        a_c_t_g_counts,size_per_cluster = get_initial_calculations(process)

    orignal_calculations_hamming_entropy = hamming_entropy_dict_per_cluster
    # temp_calculations_hamming_entropy = hamming_entropy_dict_per_cluster
    orignal_clusters = cluster_sequence_group
    # temp_cluster_sequence = cluster_sequence_group
    orignal_a_c_t_g_counts = a_c_t_g_counts
    # temp_a_c_t_g_counts = a_c_t_g_counts
    orig_size_per_cluster = size_per_cluster
    # temp_size_per_cluster = size_per_cluster


    """swapping Logic"""
    # temp_cluster_sequence = cluster_sequence_group
    iteration_info_dict = {}
    iteration_count = 1
    threshold = 800
    iterations_without_change = 0
    seq_len = len(cluster_sequence_group[list(cluster_sequence_group.keys())[0]][0][1])
    snapshot_duration = 1500
    snapshot_counter = 0

    start_time = time.process_time()
    futuretime = datetime.now()+timedelta(hours=24, seconds=1)
    # pbar = tqdm(total=1500)
    while(datetime.now()<futuretime):
    # while(iterations_without_change < threshold):
        
        temp_cluster_sequence = deepcopy(cluster_sequence_group)
        temp_calculations_hamming_entropy =  json.loads(json.dumps(hamming_entropy_dict_per_cluster))
        temp_a_c_t_g_counts = json.loads(json.dumps(a_c_t_g_counts))

        first_r_cluster,second_r_cluster,random_seq_to_swap = pick_random_seq_and_cluster(temp_cluster_sequence)

        sequence = temp_cluster_sequence[first_r_cluster][random_seq_to_swap][1]
        sequence_name = temp_cluster_sequence[first_r_cluster][random_seq_to_swap][0]

        temp_cluster_sequence[second_r_cluster].append(temp_cluster_sequence[first_r_cluster].pop(random_seq_to_swap))
        temp_size_per_cluster = {x:len(temp_cluster_sequence[x]) for x in temp_cluster_sequence.keys()}
        size_per_cluster = {x:len(cluster_sequence_group[x]) for x in cluster_sequence_group.keys()}

        """ recounting after swapping """
        temp_a_c_t_g_counts = recount_actg_after_removing(temp_a_c_t_g_counts,sequence,first_r_cluster)
        temp_a_c_t_g_counts = recount_actg_after_adding(temp_a_c_t_g_counts,sequence,second_r_cluster)
        temp_a_c_t_g_counts = calc_swap_hamming_dist_and_entropy_per_column(seq_len,temp_a_c_t_g_counts, first_r_cluster, second_r_cluster,temp_size_per_cluster)
        temp_calculations_hamming_entropy = cal_swap_hamming_entropy_per_cluster(seq_len,temp_a_c_t_g_counts,temp_size_per_cluster)
        new_en, old_en = compare_entropy_on_swap(size_per_cluster ,temp_size_per_cluster, temp_calculations_hamming_entropy,hamming_entropy_dict_per_cluster)



        if (old_en-new_en)/old_en >0.00001:
            iterations_without_change = 0
            hamming_entropy_dict_per_cluster = json.loads(json.dumps(temp_calculations_hamming_entropy))
            cluster_sequence_group = deepcopy(temp_cluster_sequence)
            a_c_t_g_counts = json.loads(json.dumps(temp_a_c_t_g_counts))
        else:
            iterations_without_change+=1
            
        
        iteration_info_dict[iteration_count]= {
                                            "seq_move":f"{first_r_cluster} --> {second_r_cluster}",
                                            "swapped_seq":sequence_name,
                                            "swap": "Failed" if (old_en-new_en)/old_en <0.00001 else "Success",
                                            "old_entropy": old_en,
                                            "new_entropy": new_en,
                                            "entropy_change": (old_en-new_en)/old_en
                                            }
        elapsed_time = timedelta(seconds =time.process_time() - start_time)
        
        # pbar.update(iteration_count)
        logger.debug(f"total_moves_made = {iteration_count}, threshold = {threshold}, total_moves_since_last_swap = {iterations_without_change}, elapsed_time: {elapsed_time},entropy({old_en, new_en})")
        # print(f"total moves made = {iteration_count} threshold = {threshold} total moves since last swap = {iterations_without_change}\nelapsed time: {elapsed_time}\nentropy({old_en, new_en})", end='\r')
        # if iteration_count % snapshot_duration == 0:
        #     print(f"total_moves_made = {iteration_count}, threshold = {threshold}, total_moves_since_last_swap = {iterations_without_change}, \nelapsed_time: {elapsed_time},entropy({old_en, new_en})")
        #     final_funtion(iteration_info_dict,
        #                     a_c_t_g_counts,orignal_a_c_t_g_counts,
        #                     orignal_clusters, cluster_sequence_group,
        #                     orignal_calculations_hamming_entropy, hamming_entropy_dict_per_cluster,
        #                     snapshot_counter = snapshot_counter
        #                     )
        #     snapshot_counter+=1
        iteration_count+=1
    logger.debug("complete")
    print(f"{process} complete")
    final_funtion(iteration_info_dict,
                            a_c_t_g_counts,orignal_a_c_t_g_counts,
                            orignal_clusters, cluster_sequence_group,
                            orignal_calculations_hamming_entropy, hamming_entropy_dict_per_cluster, process
                            )

    # return iteration_info_dict

if __name__=="__main__":
    starttime = time.time() 
    pool = mp.Pool()
    pool.map(calculate_entropy_on_swap ,range(5))
    pool.close()
    print('That took {} seconds'.format(time.time() - starttime))
    

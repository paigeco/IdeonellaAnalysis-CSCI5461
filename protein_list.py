import pandas as pd
import os
import numpy as np

#os.chdir("PET")

print("Number of Bacteria", len(os.listdir("PET")))

os.chdir("activated_deactivated")
rnaSeq_results = [fs for fs in os.listdir() if ".csv" in fs]


dict_pandas = {fs:pd.read_csv(fs) for fs in rnaSeq_results}

proteins_arr_dict = {key[:-4]:(dict_pandas[key])[['ProteinID','Annotation']].to_numpy() for key in dict_pandas.keys()}

unique_proteins = np.unique(np.array([proteins_arr_dict[key][:,0] for key in proteins_arr_dict.keys()])).tolist()
os.chdir("blastp_results")

protein_data_map = {protein:pd.read_csv(protein+".csv") for protein in unique_proteins}
"""
protein_species_list = {protein:np.unique(protein_data_map[protein][['Scientific Name']].to_numpy()) for protein in unique_proteins}
for p in protein_species_list.keys():
    if protein_species_list[p].shape[0]:
        l = protein_species_list[p]
        print(l, [x.split(" =")[0] for x in l])
"""

protein_description_count = {protein:np.array(np.unique([x.split(" [")[0].replace("TPA: ","").replace("RecName: Full=","").split(";")[0] for x in protein_data_map[protein][['Description']].to_numpy()[:,0].tolist()], return_counts=True)).T for protein in unique_proteins}
#print(protein_species_list)
for p in protein_description_count.keys():
    if protein_description_count[p].shape[0]:
        protein_description_count[p] = protein_description_count[p][np.argsort(protein_description_count[p][:,1].astype(int))[::-1]]

protein_description_most={}
for p in protein_description_count.keys():
    if protein_description_count[p].shape[0]:
        protein_description_most[p] = protein_description_count[p][:3]

print(protein_description_most)
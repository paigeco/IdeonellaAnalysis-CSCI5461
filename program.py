import pandas as pd
import numpy as np
import statistics
from math import sqrt
import scipy.stats
import seaborn
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from bs4 import BeautifulSoup
import requests



df = pd.read_csv("PET/1547922/GSE205843_PET-TPA_2Na-EG_ExpressionProfiles.txt",delimiter="\t")



# NORMALIZER
count = df[['PET_RPKM','TPA?2Na_RPKM','EG_RPKM']]

normalized_count = np.zeros(count.to_numpy().shape)

# Convert the genes column to a list
genes_in_count = count[count.columns[0]].to_list()

# Convert the list to a set (no duplicates) and count them
number_of_unique_genes = len(set(genes_in_count))

total_sample_gene_counts = [sum(count[sample].to_list()) for sample in count.columns]

median_gene_count = statistics.median(total_sample_gene_counts)

for j, row in count.iterrows():
        for i, column_name in enumerate(count.columns):
            normalized_count[j,i] = (count.loc[j,column_name]/total_sample_gene_counts[i])*median_gene_count

df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']] = normalized_count


df.to_csv("PET/1547922/GSE205843_PET-TPA_2Na-EG_ExpressionProfiles_WNORMALS.csv", index=False)

df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']] = ((df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']].to_numpy().T-np.mean(df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']].to_numpy(), axis=1).T)).T
#df[['PET_RPKM_Z','TPA?2Na_RPKM_Z','EG_RPKM_Z']] = ((df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy().T-np.mean(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=1).T)/(np.std(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=1).T+0.0000000001)).T

pet = df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy()[np.argsort(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=0)[:,0]]
tpa = df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy()[np.argsort(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=0)[:,1]]
eg = df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy()[np.argsort(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=0)[:,2]]
#print(df[['PET_RPKM_DIFF','TPA?2Na_RPKM_DIFF','EG_RPKM_DIFF']])
#print(df[['PET_RPKM_Z','TPA?2Na_RPKM_Z','EG_RPKM_Z']].to_numpy()[np.argsort(df[['PET_RPKM_Z','TPA?2Na_RPKM_Z','EG_RPKM_Z']].to_numpy(), axis=0)[:,0]])
#scipy.stats.norm.sf(

pca = PCA(n_components=2)
pca = pca.fit(df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']].to_numpy().T)

#reduced_pet = pca.transform(pet)
#reduced_tpa = pca.transform(tpa)
reduced_eg = pca.transform(df[['EG_RPKM_NORM']].to_numpy().T)
reduced_pet = pca.transform(df[['PET_RPKM_NORM']].to_numpy().T)
reduced_tpa = pca.transform(df[['TPA?2Na_RPKM_NORM']].to_numpy().T)
mean_value = pca.transform(np.mean(df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']].to_numpy(), axis=1).reshape((df.to_numpy().shape[0],1)).T)
"""
plt.title("PCA Reduced results")
plt.arrow(mean_value[0,0], mean_value[0,1],reduced_eg[0,0], reduced_eg[0,1], color="green",alpha=0.4, length_includes_head=True, head_width=500, head_length=1000)
plt.arrow(mean_value[0,0], mean_value[0,1],reduced_tpa[0,0], reduced_tpa[0,1], color="magenta",alpha=0.4, length_includes_head=True, head_width=500, head_length=1000)
plt.arrow(mean_value[0,0], mean_value[0,1],reduced_pet[0,0], reduced_pet[0,1], color="blue",alpha=0.4, length_includes_head=True, head_width=500, head_length=1000)

plt.scatter(mean_value[0,0], mean_value[0,1],s=50, color="yellow", alpha=1, label="Expression Mean")
plt.scatter(reduced_eg[0,0], reduced_eg[0,1],s=50, color="green", alpha=1, label="Ethylene Glycol")
plt.scatter(reduced_tpa[0,0], reduced_tpa[0,1],s=50, color="magenta", alpha=1, label="Terephthalic Acid")
plt.scatter(reduced_pet[0,0], reduced_pet[0,1],s=50, color="blue", alpha=1, label="Polyethylene Terephthalate")

plt.ylabel("Principal Component 1")
plt.xlabel("Principal Component 2")
#plt.xticks([])
plt.legend(loc="upper right",fontsize = 'large')
plt.show()
"""


'''
plt.title("Polyethylene Terephthalate")
plt.scatter(np.arange(0,len(pet[:,0])), pet[:,2],s=2, color="green", alpha=0.2, label="Ethylene Glycol")
plt.scatter(np.arange(0,len(pet[:,0])), pet[:,1],s=2, color="magenta", alpha=0.2, label="Terephthalic Acid")
plt.scatter(np.arange(0,len(pet[:,0])), pet[:,0],s=2, color="blue", alpha=0.9, label="Polyethylene Terephthalate")
plt.ylabel("Difference Between Counts and Average")
plt.xlabel("Sorted: By Counts Low to High")
plt.yscale('symlog')
plt.xticks([])
plt.legend(loc="upper center",fontsize = 'large')
plt.show()

plt.title("Terephthalic Acid")
plt.scatter(np.arange(0,len(tpa[:,0])), tpa[:,2],s=2, color="green", alpha=0.2, label="Ethylene Glycol")
plt.scatter(np.arange(0,len(tpa[:,0])), tpa[:,1],s=2, color="magenta", alpha=0.9, label="Terephthalic Acid")
plt.scatter(np.arange(0,len(tpa[:,0])), tpa[:,0],s=2, color="blue", alpha=0.2, label="Polyethylene Terephthalate")
plt.ylabel("Difference Between Counts and Average")
plt.xlabel("Sorted: By Counts Low to High")
plt.yscale('symlog')
plt.xticks([])
plt.legend(loc="upper center",fontsize = 'large')
plt.show()

plt.title("Ethylene Glycol")
plt.scatter(np.arange(0,len(eg[:,0])), eg[:,2],s=2, color="green", alpha=0.9, label="Ethylene Glycol")
plt.scatter(np.arange(0,len(eg[:,0])), eg[:,1],s=2, color="magenta", alpha=0.2, label="Terephthalic Acid")
plt.scatter(np.arange(0,len(eg[:,0])), eg[:,0],s=2, color="blue", alpha=0.2, label="Polyethylene Terephthalate")
plt.ylabel("Difference Between Counts and Average")
plt.xlabel("Sorted: By Counts Low to High")
plt.yscale('symlog')
plt.xticks([])
plt.legend(loc="upper center",fontsize = 'large')
plt.show()
'''
gene_expression = df[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']]

gene_expression_array = gene_expression.to_numpy()

sample_means = np.mean(gene_expression_array,axis=0)

pearson_distance = np.zeros((gene_expression_array.shape[1],gene_expression_array.shape[1]))






def get_gen_bank(ID):
    URL = """https://www.ncbi.nlm.nih.gov/protein/?term=txid1547922%5Borganism%3Aexp%5D+locus_tag%3DISF6_"""
    page = requests.get(URL+str(ID).zfill(4))
    soup = BeautifulSoup(page.content, "html.parser")
    try:
        a = (soup.find_all(class_="itemid")[0].text).split(": ")[1]
        return a
    except IndexError:
        print(str(ID)+" does not have a genome!")
        return None



def activated_deactivated(n, csv_base):
    tests = np.argsort(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=0)[:,n]
    acc =  df.take((tests[::-1])[:10])
    deac = df.take((tests[::-1])[-10:])
    acc['ProteinID'] = [get_gen_bank(x) for x in acc['ISF6_XXXX']]
    deac['ProteinID'] = [get_gen_bank(x) for x in deac['ISF6_XXXX']]
    acc.to_csv("activated_deactivated/activated_" + csv_base + ".csv", index=False)
    deac.to_csv("activated_deactivated/deactivated_" + csv_base + ".csv", index=False)
    #each = df[['ISF6_XXXX']].to_numpy()[t]
    #print("ISF6_RS"+(str(int(each*5)).zfill(5)))
    #print(prot_df[prot_df['Locus tag']=="ISF6_RS"+(str(int(each*5)).zf ill(5))],df[['Annotation']].to_numpy()[t])

def create_heat_map(n, csv_base):
    tests = np.argsort(df[['PET_RPKM_REL','TPA?2Na_RPKM_REL','EG_RPKM_REL']].to_numpy(), axis=0)[:,n]
    acc =  df.take((tests[::-1])[:10])
    deac = df.take((tests[::-1])[-10:])
    acc['ProteinID'] = [get_gen_bank(x) for x in acc['ISF6_XXXX']]
    deac['ProteinID'] = [get_gen_bank(x) for x in deac['ISF6_XXXX']]
    acc.to_csv("activated_deactivated/activated_" + csv_base + ".csv", index=False)
    deac.to_csv("activated_deactivated/deactivated_" + csv_base + ".csv", index=False)
    #each = df[['ISF6_XXXX']].to_numpy()[t]
    #print("ISF6_RS"+(str(int(each*5)).zfill(5)))
    #print(prot_df[prot_df['Locus tag']=="ISF6_RS"+(str(int(each*5)).zf ill(5))],df[['Annotation']].to_numpy()[t])


"""
activated_deactivated(0,"PolyethyleneTerephtalate")
activated_deactivated(1,"TerephtalicAcid")
activated_deactivated(2,"EthyleneGlycol")
"""
#print(get_gen_bank(528))

materials = ["PolyethyleneTerephtalate","TerephtalicAcid","EthyleneGlycol"]
def rearrange(df):
    df_out = {}
    for mk in df.keys():
        df[mk] = df[mk][['Annotation','ProteinID']]
        df_out[mk] = [row[0]+" ["+row[1]+"]" for i, row in enumerate(df[mk].to_numpy())]
    return df_out

dfs = {}
for m in materials:
    dfs[m+"_activated"] = pd.read_csv("activated_deactivated/activated_"+str(m)+".csv")
    dfs[m+"_deactivated"] = pd.read_csv("activated_deactivated/deactivated_"+str(m)+".csv")

df_for_seaborn = pd.concat(list(dfs.values()),axis=0,ignore_index=True).drop_duplicates()

import matplotlib.pyplot as plt
seaborn.set()
print(df_for_seaborn[['Annotation']].to_numpy().T[0])
seaborn.heatmap(df_for_seaborn[['PET_RPKM_NORM','TPA?2Na_RPKM_NORM','EG_RPKM_NORM']],
                yticklabels=[x[0]+" | "+x[1] for x in df_for_seaborn[['Annotation','ProteinID']].to_numpy()],
                xticklabels=["PET","TPA","EG"])
plt.show()
"""
blasting_proteins = []
for mk in dfs.keys():
    blasting_proteins += (dfs[mk]['ProteinID'].tolist())

print(set(blasting_proteins))
"""

for x, column_name_x in enumerate(gene_expression.columns):
    for y, column_name_y in enumerate(gene_expression.columns):
        if x == y:
            pearson_distance[x,y] = float(0)
        else:
            x_bar = gene_expression_array[:,x] - sample_means[x]
            y_bar = gene_expression_array[:,y] - sample_means[y]
            pearson_distance[x,y] = np.sum(np.multiply(x_bar,y_bar)) / (sqrt(np.sum(np.square(x_bar)))*sqrt(np.sum(np.square(y_bar))))

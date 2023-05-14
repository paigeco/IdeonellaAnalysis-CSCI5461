import pandas as pd
import numpy as np
import os
from bs4 import BeautifulSoup
import requests

def download_proteins_list(URL, k, plastic_type):
    page = requests.get(URL)
    soup = BeautifulSoup(page.content, "html.parser")
    try:
        a = soup.find_all(class_="GenomeList2")[0].find_all("tr",align="center")
    except IndexError:
        print(str(k)+" does not have a genome!")
        return
    b = list(a[0].find_all("td")[6].children)[0].get_attribute_list("href")[0]
    c = str(b).split("/proteins/")[1].split("|")[0].split("/")
    
    ids = ""
    for x in a:
        y = x.find_all("td")
        c_name = y[0].find("span").get_attribute_list("title")[0] 
        c_num = y[1].get_text()
        ids = ids + '''"''' + c_name + ("" if c_num == "-" else " "+str(c_num)) + '''",'''
    
    payload = {
        "q":('''[display()].from(GenomeProteins).usingschema(/schema/GenomeAssemblies).matching(genome_id==["'''+str(c[0])+'''"] and genome_assembly_id==["'''+str(c[1])+'''"] and replicon_name==[''' + ids[:-1] + '''])'''),
        "fields":"replicon_name|Name,replicon_accession|Accession,start|Start,stop|Stop,strand|Strand,gene_id|GeneID,locus|Locus,locus_tag|Locus tag,accession|Protein product,length|Length,name|Protein Name",
        'filename':'proteins_'+str(c[0])+'_'+str(c[1])+'.csv','nolimit':'on'
    }
    url = "https://www.ncbi.nlm.nih.gov/genomes/solr2txt.cgi"
    r = requests.get(url, params=payload, allow_redirects=True)
    open(plastic_type+'/'+str(k)+'/proteins.csv', 'wb').write(r.content)
    #quit()


df = pd.read_csv("polymers.tsv",delimiter="\t")


data = df.to_numpy()[:,:2]

data[:,0] = [x[:-1] for x in data[:,0]]

polymer_name_map = {each[0]:each[1] for each in data}

# read the list of degraders
df = pd.read_csv("degraders_list.tsv",delimiter="\t")

# get list of unique plastics
unique_plastics = np.unique(df['Plastic'].to_numpy())


plastics_to_include = ["PET","PETG"]
plastics_to_exclude = ["HDPE","LDPE","LDPE Blend","LLDPE","LLDPE Blend","Nylon","PC","PE","PE Blend","PS","PS Blend","PVC","PVC Blend"]

is_PET = np.sum(np.array([(df['Plastic'] == plastic).to_numpy() for plastic in plastics_to_include]),axis=0).astype(bool)
is_OTH = np.sum(np.array([(df['Plastic'] == plastic).to_numpy() for plastic in plastics_to_exclude]),axis=0).astype(bool)

a = df[is_PET]
b = a[np.array([not ( "sp." in x ) for x in a['Microorganism'].to_numpy()])]
c = b[np.array([( "superkingdom:Bacteria" in x ) for x in b['lineage'].to_numpy()])]
pet_df = c

a = df[is_OTH]
b = a[np.array([not ( "sp." in x ) for x in a['Microorganism'].to_numpy()])]
c = b[np.array([( "superkingdom:Bacteria" in x ) for x in b['lineage'].to_numpy()])]
other_df = c

#(df[is_OTH]).to_csv("OTHDegs")
#print(np.unique(c['Tax ID']))
#degrades_pet = pet_df[['Tax ID','Microorganism','Enzyme','Plastic','Enzyme']].drop_duplicates().sort_values(by=['Tax ID'])
#degrades_pyrolysis = other_df[['Tax ID','Microorganism','Enzyme','Plastic']].drop_duplicates().sort_values(by=['Tax ID'])


degrade = {}
for x in np.unique(pet_df.to_numpy()[:,0]):
    degrade[x] = []


for k in np.unique(pet_df['Tax ID'].to_numpy()):
    try:
        os.mkdir('PET/'+str(k))
    except FileExistsError:
        pass
    print(k)
    (pet_df[pet_df['Tax ID']==k]).to_csv("PET/"+str(k)+"/info.csv",index=False)
    download_proteins_list("https://www.ncbi.nlm.nih.gov/genome/?term=txid"+str(k)+"[Organism:exp]",k,"PET")




for k in np.unique(other_df['Tax ID'].to_numpy()):
    try:
        os.mkdir('OTHER/'+str(k))
    except FileExistsError:
        pass
    print(k)
    (other_df[other_df['Tax ID']==k]).to_csv("OTHER/"+str(k)+"/info.csv",index=False)
    download_proteins_list("https://www.ncbi.nlm.nih.gov/genome/?term=txid"+str(k)+"[Organism:exp]",k,"OTHER")

#"document.getElementsByClassName("GenomeList2")[0].getElementsByTagName("tr")[1].getElementsByTagName("td")[6].firstChild.href"


# degrades_pet.to_csv('pet.csv', index=False)
# degrades_pyrolysis.to_csv('pyr.csv', index=False)
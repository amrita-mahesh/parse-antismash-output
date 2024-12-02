"""
Parse AntiSMASH output files

Note: This code is compatible with output .json files produced by AntiSMASH 7.0.

The output excel file has 2 sheets: Sheet 1 tabulates the genome name, cluster
type (e.g.: terpene, RiPP-like) and count of each cluster type. Sheet 2 tabulates all
the identified BGCs which have a match in the MiBiG database.

"""

import json
import os
import pandas as pd


#Change working directory to directory with all the .json files produced by AntiSMASH.
working_dir = ""

os.chdir(working_dir)
json_files = [x for x in os.listdir() if x.endswith('json')]

antismash_data = []
known_clusters_data = []

for input_file in json_files:
    genome_name = input_file.replace('.json', '')
   
    f = open(input_file)
    data = json.load(f)
    f.close()
    
    #Identify records with non-zero number of clusters
    records = data['records']
    valid_records = []
    cluster_count = {}
    
    for record in records:
        regions = record['areas']
        if len(regions) != 0:
            valid_records.append(record)         
            
    #Find the set of clusters/known clusters for all records  
    
    for record in valid_records:  
        regions = record['areas']
        known_clusters = record['modules']['antismash.modules.clusterblast']['knowncluster']['results']
        
        assert(len(regions) == len(known_clusters))
    
        for region in regions:
            product = ', '.join(region['products'])
            if product in cluster_count.keys():
                cluster_count[product]+=1
            else:
                cluster_count[product] = 1
                
        #Extract all clusters with a match in MiBiG
        for cluster in known_clusters:
            if cluster['total_hits'] != 0:
                known_cluster_match = cluster['ranking'][0]
                match_name = known_cluster_match[0]['description']
                match_type = known_cluster_match[0]['cluster_type']
                match_similarity = known_cluster_match[1]['similarity']
                
                known_clusters_data.append([genome_name, product, match_name, match_type, match_similarity])
        
    for x in cluster_count:
        antismash_data.append([genome_name, x, cluster_count[x]])


#Output data to excel files

output_1 = pd.DataFrame(antismash_data, columns = ['Genome', 'Gene_cluster', "Count"])
output_2 = pd.DataFrame(known_clusters_data, columns = ['Genome', 'Gene_cluster','Known_Cluster', "Cluster_type", 'Similarity'])

with pd.ExcelWriter('antismash_outputs.xlsx') as writer:
    output_1.to_excel(writer, sheet_name='Gene_clusters', index = False)
    output_2.to_excel(writer, sheet_name='Known_clusters', index = False)

     
            
    
    
            


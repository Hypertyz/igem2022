# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:20:56 2022

@author: HoangPham
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import re
import pandas as pd
# writing to file
# Using readlines()
file1 = open('GCA_000092025.1_ASM9202v1_cds_from_genomic.fna', 'r')
Lines = file1.readlines()


protein_id = []
protein = []
chromosome = []
location = []
gene = []
gene_id = []
# Strips the newline character
for idx, line in enumerate(Lines):
    if ">" in line:
        chromosome.append(line[5:15])
        desc = re.findall('\[(.*?)\]', line)
        if len(desc) == 5:
                gene_name = "no name"
                gene.append(gene_name)
        if len(desc) == 6:
            gene_name = desc[0][5:]
            gene.append(gene_name)
        for item in desc:
            if "protein_id=" in item:
                pro_id = item[11:]
                protein_id.append(pro_id)
            if "protein=" in item:
                pro = item[8:]
                protein.append(pro)
            if "location=" in item:
                loc = item[9:]
                location.append(loc)
            if "locus_tag=" in item:
                locus = item[10:]
                gene_id.append(locus)
            
            
metadata = pd.DataFrame(list(zip(protein_id, gene_id, gene, protein, chromosome, location)), 
                        columns = ['Protein_id', 'Gene_id', 'Gene_name', 'Protein', 'Chromosome', 'Location'])

metadata.to_csv('reads_metadata.csv', index=False)

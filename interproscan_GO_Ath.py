#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :zhangyan1220zj@163.com
# @FileName :interproscan_GO_Ath.py
# @Time :2023/11/22 16:07
# @Last time: 2025/6/25 22:49
# python3 interproscan_GO_Ath.py
# Get the GO number and the associated GO function based on the interproscan results.
import os
import re
import glob
import time
import shutil
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


# # 1 Generate files with GO numbers
def process_file(input_file, output_folder):
    with open(input_file, 'r') as infile:
        file_name = os.path.basename(input_file)
        output_file_path = os.path.join(output_folder, file_name)
        with open(output_file_path, 'w') as outfile:
            for line in infile:
                if "GO:" in line:
                    new_line = line.split("\t")
                    id, length, database, database_hao, Description, Start, End, Score, InterPro_Accession, InterPro_Description, GO_Terms = new_line[0], new_line[2], new_line[3], new_line[4], new_line[5], new_line[6], new_line[7], new_line[8], new_line[11], new_line[12], new_line[13]
                    outfile.write(f"{id}\t{length}\t{database}\t{database_hao}\t{Description}\t{Start}\t{End}\t{Score}\t{InterPro_Accession}\t{InterPro_Description}\t{GO_Terms}\n")


# 2 Correspondence between gene id and GO number
def process_go_file(input_file, output_folder):
    data_dict = {}
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            key, e_value, value = parts[0], parts[7], parts[-1]
            if key in data_dict and e_value < data_dict[key][1]:
                data_dict[key] = [value, e_value]
            elif key not in data_dict:
                data_dict[key] = [value, e_value]

    file_name = os.path.basename(input_file)
    new_file_path = os.path.join(output_folder, file_name)
    with open(new_file_path, 'w') as new_file:
        for key, value in data_dict.items():
            new_file.write(f"{key}\t{value[0]}\n")


# 3 The result file from the second step was further processed into, GO-gene id, number of genes, GO name, type it belongs to.
def create_go_terms_dict(file_name):
    dict_go = {}
    with open(file_name, 'r') as go_file:
        for line in go_file:
            parts = line.split("\t")
            if len(parts) >= 3:
                dict_go[parts[0]] = parts[1:3]
    return dict_go
def process_output_go_files(input_folder, output_folder, go_terms_file):
    go_anno = create_go_terms_dict(go_terms_file)
    file_list = glob.glob(os.path.join(input_folder, "*"))
    for file_path in file_list:
        go_dict = defaultdict(list)
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                gene_id = parts[0]
                go_terms = parts[1].split('|')
                for go_term in go_terms:
                    go_dict[go_term].append(gene_id)
        output_file_path = os.path.join(output_folder, os.path.basename(file_path))
        with open(output_file_path, 'w') as output_file:
            for go_term, gene_ids in go_dict.items():
                seq = ','.join(gene_ids)
                go_type = "\t".join(go_anno.get(go_term, ["", ""]))
                output_file.write(f"{go_term}\t{seq}\t{len(gene_ids)}\t{go_type}\n")


# 4 Pairs of event-related genes, changed to one column and de-weighted
def process_single_file(file_path, output_folder):
    gene_info_list = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            gene_info_list.append(f"{parts[1]}")
            gene_info_list.append(f"{parts[0]}")

    unique_gene_info = list(set(gene_info_list))
    new_name = os.path.basename(file_path)
    output_file_path = os.path.join(output_folder, new_name)
    with open(output_file_path, 'w') as output_file:
        output_file.write('\n'.join(unique_gene_info))


# 5 GO_Terms for event-related genes
def process_event_genes_with_GO(output_4_event_genes_folder, output_2_gene_GO_folder, output_5_event_genes_GO):
    shutil.rmtree(output_5_event_genes_GO, ignore_errors=True)
    os.makedirs(output_5_event_genes_GO)

    for file_path in glob.glob(os.path.join(output_4_event_genes_folder, "*")):
        parts = os.path.basename(file_path).split("_")[0]
        corresponding_file_path = os.path.join(output_2_gene_GO_folder, f"{parts}.pep.tsv")

        gene_info_dict = {}
        with open(corresponding_file_path, 'r') as corresponding_file:
            gene_info_dict = {key: value for key, value in (line.strip().split('\t')[:2] for line in corresponding_file)}

        with open(file_path, 'r') as file_name_open, open(f"{output_5_event_genes_GO}/{os.path.basename(file_path)}", 'w') as new_file:
            for line in file_name_open:
                new_line = line.strip()
                if new_line in gene_info_dict:
                    new_file.write(f"{new_line}\t{gene_info_dict[new_line]}\n")


# 6 The result file from step 5, is further processed into, GO-gene id, number of genes, GO name, type it belongs to. The call is to the code from step 5
# 7 GO enrichment analysis, numerical values, P-value, etc.
def process_event_genes_with_GO_pvalue(output_6_event_GO_gene, output_3_GO_gene, output_7_event_genes_GO_gene):
    shutil.rmtree(output_7_event_genes_GO_gene, ignore_errors=True)
    os.makedirs(output_7_event_genes_GO_gene)

    lens_files = glob.glob(f"input_lens/*")
    all_gene_counts = {os.path.basename(file).split(".")[0]: pd.read_csv(file, sep='\t', header=None).iloc[:, 1].sum() for file in lens_files}

    event_files = glob.glob(f"output_4_event_genes/*")
    event_gene_counts = {os.path.basename(file).split(".")[0]: pd.read_csv(file, sep='\t', header=None).shape[0] for file in event_files}
    # print(event_gene_counts, output_6_event_GO_gene, output_3_GO_gene, output_7_event_genes_GO_gene)
    for file_path in glob.glob(os.path.join(output_6_event_GO_gene, "*")):
        new_part = os.path.basename(file_path)
        parts = new_part.split("_")[0]
        corresponding_file_path = os.path.join(output_3_GO_gene, f"{parts}.pep.tsv")
        # print(corresponding_file_path)

        gene_info_dict = {}
        with open(corresponding_file_path, 'r') as corresponding_file:
            # gene_info_dict = {key: value for key, _, value,_ , _ in (line.strip().split('\t') for line in corresponding_file)}
            # print(corresponding_file)
            for line in corresponding_file:
                columns = line.strip().split('\t')
                key = columns[0]
                value = columns[2]
                gene_info_dict[key] = value

        with open(file_path, 'r') as file_name_open, open(f"{output_7_event_genes_GO_gene}/{os.path.basename(file_path)}", 'w') as new_file:
            for line in file_name_open:
                new_line_hh = line.strip().split("\t")[0]
                seq = "\t".join(line.strip().split("\t")[1:3])
                sample_number = line.strip().split("\t")[2]
                other = "\t".join(line.strip().split("\t")[3:])

                total_genes = all_gene_counts[parts]  # Total number of genes
                sample_genes = event_gene_counts[new_part.split(".")[0]]  # Total number of genes in the sample

                if new_line_hh in gene_info_dict.keys():
                    known_go_genes = gene_info_dict[new_line_hh]  # The number of genes in the GO term is known
                    significant_genes = sample_number  # Number of genes significantly enriched in the sample
                    p_value = hypergeom.sf(int(significant_genes) - 1, int(total_genes), int(known_go_genes), int(sample_genes))
                    new_file.write(f"{new_line_hh}\t{seq}\t{sample_genes}\t{known_go_genes}\t{total_genes}\t{p_value}\t{other}\n")


# 8.1 Complete Benjamini-Hochberg Correction Functions
def bh_correction(p_values):
    m = len(p_values)
    sorted_idx = p_values.argsort()
    reverse_idx = np.empty_like(sorted_idx)
    reverse_idx[sorted_idx] = np.arange(m)
    sorted_p = p_values[sorted_idx]
    q_values = sorted_p * m / (np.arange(m) + 1)
    q_min = q_values[::-1].cummin()[::-1]
    return np.minimum(q_min, 1.0)[reverse_idx]

# 8.2 Perform p-value correction on individual files and save the results
def adjust_p_values_in_file(input_path, output_path):
    data = pd.read_csv(input_path, sep='\t', header=None)
    p_values = data[6]  # Column 7 is the p-value
    q_values = bh_correction(p_values)
    columns = data.columns.tolist()
    new_columns = columns[:7] + ["q_value"] + columns[7:]
    data.insert(7, "q_value", q_values)
    data.columns = new_columns
    data.to_csv(output_path, sep='\t', index=False, header=False)

# 8.3 Processing of files in an entire folder
def process_folder_step8(input_folder, output_folder):
    input_files = glob.glob(os.path.join(input_folder, '*'))
    print(input_files)
    for input_file in input_files:
        if os.path.isfile(input_file):  
            filename = os.path.basename(input_file).split(".")[0] + "_adjust_pvalue.txt"
            output_path = os.path.join(output_folder, filename)
            # try:
            adjust_p_values_in_file(input_file, output_path)
            print(f"Step 8 over: {filename} -> {output_folder}")
            # except Exception as e:
            #     print(f"Step 8 error: {filename}: {str(e)}")


# 9.1 Convert GO enrichment analysis results to R-compatible formats
def convert_to_r_format(input_file_path, output_file_path):
    # try:
    data = pd.read_csv(input_file_path, sep='\t', header=None)
    if len(data.columns) < 10:
        print(f"Warning: file {os.path.basename(input_file_path)} needs at least 10 columns of data, skip conversion")
        return False
    converted_data = pd.DataFrame()
    ontology_mapping = {
        'molecular_function': 'MF',
        'biological_process': 'BP',
        'cellular_component': 'CC'
    }
    print(data[9])
    # converted_data['ONTOLOGY'] = data[9].apply(
    #     lambda x: ontology_mapping.get(x.strip(), 'Unknown')
    # )
    converted_data['ONTOLOGY'] = data[9].apply(
        lambda x: ontology_mapping.get(str(x).strip(), 'Unknown') 
        if pd.notna(x) else ontology_mapping.get('NaN', 'Unknown')
    )
    converted_data['ID'] = data[0]
    converted_data['Description'] = data[8]
    converted_data['GeneRatio'] = data[2].astype(str) + '/' + data[3].astype(str)
    converted_data['BgRatio'] = data[4].astype(str) + '/' + data[5].astype(str)
    converted_data['pvalue'] = data[6]
    converted_data['p.adjust'] = data[6]
    converted_data['qvalue'] = data[7]
    converted_data['geneID'] = data[1]
    converted_data['Count'] = data[2]
    converted_data.to_csv(output_file_path, sep='\t', index=False)
    return True
    # except Exception as e:
    #     print(f"Error while converting file {os.path.basename(input_file_path)}: {str(e)}")
    #     return False

# 9.2 Processing of GO enrichment analysis files from entire folders
def process_folder_step9(input_folder, output_folder):
    input_files = glob.glob(os.path.join(input_folder, '*'))
    processed_count = 0
    skipped_count = 0
    for input_file_path in input_files:
        if os.path.isfile(input_file_path):
            filename = os.path.basename(input_file_path).split(".")[0] + "_runR.txt"
            output_file_path = os.path.join(output_folder, filename)

            if convert_to_r_format(input_file_path, output_file_path):
                processed_count += 1
                print(f"✓ over: {filename}")
            else:
                skipped_count += 1


if __name__ == "__main__":
    input_interproscan_tsv = "input_interproscan_tsv"  # (changeable) - input - interprosca results file folder
    input_event_genes = "input_event_genes"  # (changeable) - Input - event related gene pair files
    output_1_haveGO_line = "output_1_haveGO_line"  # (changeable) - Output - folder with Go message lines
    output_2_gene_GO = "output_2_gene_GO"  # (changeable) - Output - folder with information about the GO number corresponding to the gene id
    output_3_GO_gene = "output_3_GO_gene"  # (changeable) - Output - GO, gene id, number of genes, GO name, type
    output_4_event_genes = "output_4_event_genes"  # (changeable) - Output - single column of de-duplicated files
    output_5_event_genes_GO = "output_5_event_genes_GO"  # (changeable) - Output - GO files for event-related genes
    output_6_event_GO_gene  = "output_6_event_GO_gene"  # (changeable) - Output - GO, gene id, number of genes, GO name, type
    output_7_event_genes_GO_gene = "output_7_event_genes_GO_gene"  # (subject to change) - Output - number counts, p-value
    output_8_adjust_pvalue = "output_8_adjust_pvalue"  # (changeable) - output corrected p-value

    go_terms_file = "go_terms.txt"  # Change - the function name file corresponding to the GO number
    start_time = time.time()

    # 1 Generate files with GO numbers
    shutil.rmtree(output_1_haveGO_line, ignore_errors=True)
    os.makedirs(output_1_haveGO_line)
    input_files = glob.glob(os.path.join(input_interproscan_tsv, "*.*"))
    with ProcessPoolExecutor(max_workers=8) as executor:  # Change - maximum number of threads
        executor.map(process_file, input_files, [output_1_haveGO_line] * len(input_files))

    # 2 Correspondence between gene id and GO number
    shutil.rmtree(output_2_gene_GO, ignore_errors=True)
    os.makedirs(output_2_gene_GO)
    file_list = glob.glob(os.path.join(output_1_haveGO_line, '*.tsv'))
    for file_path in file_list:
        process_go_file(file_path, output_2_gene_GO)
    
    # 3 The result file from the second step was further processed into, GO-gene id, number of genes, GO name, type it belongs to.
    shutil.rmtree(output_3_GO_gene, ignore_errors=True)
    os.makedirs(output_3_GO_gene)
    process_output_go_files(output_2_gene_GO, output_3_GO_gene, go_terms_file)

    # 4 Pairs of event-related genes, changed to one column and de-weighted
    shutil.rmtree(output_4_event_genes, ignore_errors=True)
    os.makedirs(output_4_event_genes)
    file_list = glob.glob(os.path.join(input_event_genes, "*"))
    for file_path in file_list:
        process_single_file(file_path, output_4_event_genes)

    # 5 GO_Terms for event-related genes
    shutil.rmtree(output_5_event_genes_GO, ignore_errors=True)
    os.makedirs(output_5_event_genes_GO)
    process_event_genes_with_GO(output_4_event_genes, output_2_gene_GO, output_5_event_genes_GO)

    # 6 The result file from step 5, was further processed into, GO-gene id, number of genes, GO name, type it belongs to.
    shutil.rmtree(output_6_event_GO_gene, ignore_errors=True)
    os.makedirs(output_6_event_GO_gene)
    process_output_go_files(output_5_event_genes_GO, output_6_event_GO_gene, go_terms_file)

    # 7 GO enrichment analysis, numerical values, P-values, etc.
    process_event_genes_with_GO_pvalue(output_6_event_GO_gene, output_3_GO_gene, output_7_event_genes_GO_gene)

    # 8 P-value correction using the Benjamini-Hochberg correction function
    shutil.rmtree(output_8_adjust_pvalue, ignore_errors=True)
    os.makedirs(output_8_adjust_pvalue)
    process_folder_step8(output_7_event_genes_GO_gene, output_8_adjust_pvalue)

    # 9.2 Processing GO enrichment analysis files throughout the folder
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total time taken: {elapsed_time:.2f} seconds")


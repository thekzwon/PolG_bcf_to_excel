#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 11:20:59 2020

@author: kendra
"""

# Final python script to convert bcfquery files into excel

# Updated 9_15_20

import csv
from os import listdir
from collections import Counter
from Bio import SeqIO

path = "/Users/kendra/Desktop/Research/Havird_Lab/PolG/total_python_input/"

####################################################################
####Codon_table#######
def translate(seq):

    table = {
        'ATA': 'M', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': '_', 'AGG': '_',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'W', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

###takes single letter amino acid abbreviation and returns group###
def hydrophobicity(amino_acid):
    table = {
        'A': 'Hydrophobic', 'R': 'Basic', 'N': 'Hydrophilic',
        'D': 'Acidic', 'C': 'Hydrophilic', 'E': 'Acidic',
        'Q': 'Hydrophilic', 'G': 'Hydrophilic', 'H': 'Basic',
        'I': 'Hydrophobic', 'L': 'Hydrophobic', 'K': 'Basic',
        'M': 'Hydrophobic', 'F': 'Hydrophobic', 'P': 'Hydrophobic',
        'S': 'Hydrophilic', '_': 'Stop', 'T': 'Hydrophilic',
        'W': 'Hydrophobic', 'Y': 'Hydrophobic', 'V': 'Hydrophobic',
    }
    return(table[amino_acid])

###Takes the complement of bases###
def complement(bases):
    complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    complement_str = ""
    for nucleotide in bases:
        complement_str += (complement_dict[nucleotide])
    return(complement_str)

##################################################################

# working point mutation list
point_mutation_list = []
# working indel list
indel_mutation_list = []
#Some mutations were printed as CA ->TA rather than C->T
#This list processes those mutations
fake_indel_list=[]
# list that gets added to file
final_mut_list = []




#List for renaming animals to lowercase letters
new_animal_names = [chr(i) for i in range(97,123)]


# Gets and prints list of files in directory
bcf_files = [file for file in listdir(path)]
print(bcf_files)

animal_num = 0
sample_num = 0

# Clealing and combining files
for file in sorted(bcf_files):
    # Can change query to whatever the files are named
    if file[-5:] == "query":
        open_file = open(path + file, "r")
        read_file = open_file.readlines()
        for line in read_file:
            line = line.replace(
                "    ",
                " ").replace(
                "   ",
                " ").replace(
                "\n",
                "").replace("\t"," ").split(" ")
                    
            #Adds the unique file information to each line of the file.
            animal_name = file[:3].strip("_")[:-1]
            tissue = file[:3].strip("_")[-1]
            sample_name = file[:3].strip("_")
            raw_mutation_line = []
            raw_mutation_line.append(animal_name)
            raw_mutation_line.append(new_animal_names[animal_num])
            raw_mutation_line.append(
                tissue.replace(
                    "B", "Brain").replace(
                    "L", "Liver"))
            raw_mutation_line.append(sample_name)
            raw_mutation_line = raw_mutation_line + line
            # finding alternate alleles, may be multiple
            alternate_alleles = raw_mutation_line[-1].split(",")
            lines_to_add = len(alternate_alleles)

            # adding lines for each alternate allele
            while lines_to_add > 0:
                mutation_line = []
                mutation_line = raw_mutation_line[0:6]
                index = lines_to_add
                allele_frequency = raw_mutation_line[6].split(",")
                count = raw_mutation_line[7].split(",")
                alt_allele = raw_mutation_line[8].split(",")
                mutation_line.append(allele_frequency[index - 1])
                mutation_line.append(count[0])
                mutation_line.append(count[index])
                mutation_line.append(alt_allele[index - 1])


                lines_to_add = lines_to_add - 1
                # separating snps and indels
                if len(mutation_line[5]) == 1 and len(mutation_line[9]) == 1:
                    point_mutation_list.append(mutation_line)
                else:
                    if len(mutation_line[5]) == len(mutation_line[9]):
                        mut_num = 0
                        for base in mutation_line[5]:
                            if mutation_line[5][mut_num] != mutation_line[9][mut_num]:
                                mut_line = []
                                mut_line += mutation_line
                                mut_line[5] = mut_line[5][mut_num]
                                mut_line[9] = mut_line[9][mut_num]
                                fake_indel_list.append(mut_line)
                            mut_num += 1
                    else:
                        indel_mutation_list.append(mutation_line)

        sample_num += 1
        if sample_num % 2 == 0:
            animal_num += 1


##Adds fake indels to total sheet. If the mutation is a duplicate, adds to previous mutation.##
            
repeat_det=[]
for mutation in fake_indel_list:
    repeat = [mutation[3],mutation[4], mutation[5],mutation[9]]
    repeat_det.append(repeat)            

indel_dups = []
repeat_dups = []
dup_count = -1

for mutation in fake_indel_list:
    repeat = [mutation[3],mutation[4], mutation[5],mutation[9]]
    if repeat_det.count(repeat)==1:
        point_mutation_list.append(mutation)
    else:
        if repeat_dups.count(repeat)==0:
            repeat_dups.append(repeat)
            indel_dups.append(mutation)
            dup_count+=1
        else:
            indel_dups[dup_count][6]=str(float(mutation[6])+float(indel_dups[dup_count][6]))
            indel_dups[dup_count][8]=int(mutation[8])+ int(indel_dups[dup_count][8])
        



som_germ_det = []
for mutation in point_mutation_list:
    som_germ = mutation[0] + mutation[4] + mutation[5] + mutation[9]
    som_germ_det.append(som_germ)
somgermcounter = (Counter(som_germ_det))


for mutation in point_mutation_list:
    som_germ = mutation[0] + mutation[4] + mutation[5] + mutation[9]
    mutation.append(str(somgermcounter[som_germ]).replace(
     "1", "Somatic").replace("2", "Germline"))

point_mutation_list[:] = (value for value in point_mutation_list if value != "")


som_germ_indel = []
for mutation in indel_mutation_list:
    som_germ = mutation[0] + mutation[4] + mutation[5] + mutation[9]
    som_germ_indel.append(som_germ)

somgermcounter = (Counter(som_germ_indel))


for mutation in indel_mutation_list:
    som_germ = mutation[0] + mutation[4] + mutation[5] + mutation[9]
    mutation.append(str(somgermcounter[som_germ]).replace(
        "1", "Somatic").replace("2", "Germline"))
    if len(mutation[5]) > len(mutation[9]):
        mutation.append("deletion")
    elif len(mutation[5]) == len(mutation[9]):
        mutation.append("NA")
    else:
        mutation.append("insertion")


for mutation in indel_mutation_list:
    if len(mutation[5]) - len(mutation[9]) % 3 == 0:
        mutation.append("no frameshift")
    else:
        mutation.append("frameshift")
        

seq_record = list(SeqIO.parse(open(path + "sequence.gb"), "gb"))

# open genbank file and parse
gene_records = []


for record in seq_record:
    for feature in record.features:
        single_gene = []
        if feature.type == "CDS":
            single_gene.append(str(feature.qualifiers['gene'])[2:-2])
            single_gene.append(feature.type)
            single_gene.append(feature.location._start.position)
            single_gene.append(feature.location._end.position)
            single_gene.append(feature.location.strand)
            single_gene.append(str(feature.location.extract(record.seq)))
            gene_records.append(single_gene)
        elif feature.type == "rRNA" or feature.type == "tRNA" or feature.type == "D-loop":
            single_gene.append(feature.type)
            single_gene.append(feature.type)
            single_gene.append(feature.location._start.position)
            single_gene.append(feature.location._end.position)
            single_gene.append(feature.location.strand)
            single_gene.append(str(feature.location.extract(record.seq)))
            gene_records.append(single_gene)


for record in seq_record:
    mtDNA_seq = record.seq

rRNA_len = 0
tRNA_len = 0
dloop_len = 0
CDS_len = 0

for gene in gene_records:
    if gene[1] == "rRNA":
        rRNA_len += gene[3] - gene[2]
    if gene[1] == "tRNA":
        tRNA_len += gene[3] - gene[2]
    if gene[1] == "D-loop":
        dloop_len += gene[3] - gene[2]
    if gene[1] == "CDS":
        CDS_len += gene[3] - gene[2]


region_length_dictionary = {
    "rRNA": rRNA_len,
    "tRNA": tRNA_len,
    "D-loop": dloop_len,
    "CDS": CDS_len}


    
    
for gene in gene_records:
    #trNA, rRNA, D-loop
    if gene[1] != "CDS":
        gene_length = gene[3] - gene[2]
        for mutation in point_mutation_list:
            final_mutation = []
            final_mutation += mutation
            mutation_loc = mutation[4]
            if int(mutation_loc) > gene[2] and int(mutation_loc) <+ gene[3]:
                ref = int(mutation_loc) - gene[2]
                gene_seq = gene[-1]
                surrounding = mtDNA_seq[int(
                    mutation[4]) - 2] + mtDNA_seq[int(mutation_loc)]
                final_mutation.append(str(gene[0]))
                final_mutation.append(str(gene[0]))
                final_mutation.append(surrounding)
                final_mutation.append(
                    str(float(final_mutation[6]) / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(float(final_mutation[6]) / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(1 / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(1 / region_length_dictionary[gene[1]]))
                final_mutation.append("NA")
                final_mutation.append("NA")
                final_mutation.append("NA")
                final_mutation.append("NA")
                final_mutation.append("NA")
                final_mutation.append("NA")
                final_mut_list.append(final_mutation)


# Forward CDS
    elif gene[-2] == 1:
        gene_length = gene[3] - gene[2]
        for mutation in point_mutation_list:
            final_mutation = []
            final_mutation += mutation
            mutation_loc = mutation[4]
            if int(mutation_loc) > gene[2] and int(mutation_loc) <= gene[3]:
                ref = int(mutation_loc) - gene[2]
                # adding polyA tail: needed for some stop codons
                gene_seq = gene[-1] + "AA"
                surrounding = str(gene_seq[ref - 2]) + str(gene_seq[ref])
                b = (ref // 3)
                # Gets the correct codon by using the reference sequence and
                # translates
                right_codon = (gene_seq[(b * 3):(b * 3) + 3])
                codon_list = [right_codon[0], right_codon[1], right_codon[2]]
                # Gets the new(mutated) codon
                codon_list[ref % 3] = mutation[5]
                new_codon = ''
                new_codon = new_codon.join(codon_list)
                original_aa = translate(right_codon.upper())
                mutated_aa = translate(new_codon.upper())
                final_mutation.append(str(gene[0]))
                final_mutation.append(str(gene[1]))
                final_mutation.append(surrounding)
                final_mutation.append(
                    str(float(final_mutation[6]) / gene_length))
                final_mutation.append(str(float(final_mutation[6]) / CDS_len))
                final_mutation.append(str(1 / gene_length))
                final_mutation.append(str(1 / CDS_len))
                final_mutation.append(original_aa)
                final_mutation.append(hydrophobicity(original_aa))
                final_mutation.append(mutated_aa)
                final_mutation.append(hydrophobicity(mutated_aa))

                if final_mutation[18] == final_mutation[20]:
                    final_mutation.append("silent")
                    final_mutation.append("NA")
                elif final_mutation[20] == "_":
                    final_mutation.append("nonsense")
                    final_mutation.append("radical")
                else:
                    final_mutation.append("missense")
                    if final_mutation[19] == final_mutation[21]:
                        final_mutation.append("conservative")
                    else:
                        final_mutation.append("radical")
                final_mut_list.append(final_mutation)


# Reverse CDS
    else:
        gene_length = gene[3] - gene[2]
        for mutation in point_mutation_list:
            final_mutation = []
            final_mutation += mutation
            mutation_loc = mutation[4]
            if int(mutation_loc) > gene[2] and int(mutation_loc) <= gene[3]:
                ref = int(mutation_loc) - gene[3]

                surrounding = complement(
                    str(gene_seq[ref - 2]) + str(gene_seq[ref]))

                # adding polyA tail: needed for some stop codons
                gene_seq = gene[-1] + "AA"
                b = (abs(ref) // 3)
                # Gets the correct codon by using the reference sequence and
                # translates
                right_codon = (gene_seq[(b * 3):(b * 3) + 3])
                codon_list = [right_codon[0], right_codon[1], right_codon[2]]

                # Gets the new(mutated) codon

                codon_list[ref % 3] = complement(mutation[5])
                new_codon = ''
                new_codon = new_codon.join(codon_list)
                original_aa = translate(right_codon.upper())
                mutated_aa = translate(new_codon.upper())
                final_mutation.append(str(gene[0]))
                final_mutation.append(str(gene[1]))
                final_mutation.append(surrounding)
                final_mutation.append(
                    str(float(final_mutation[6]) / gene_length))
                final_mutation.append(str(float(final_mutation[6]) / CDS_len))
                final_mutation.append(str(1 / gene_length))
                final_mutation.append(str(1 / CDS_len))
                final_mutation.append(original_aa)
                final_mutation.append(hydrophobicity(original_aa))
                final_mutation.append(mutated_aa)
                final_mutation.append(hydrophobicity(mutated_aa))

                if final_mutation[18] == final_mutation[20]:
                    final_mutation.append("silent")
                    final_mutation.append("NA")
                elif final_mutation[20] == "_":
                    final_mutation.append("nonsense")
                    final_mutation.append("radical")
                else:
                    final_mutation.append("missense")
                    if final_mutation[19] == final_mutation[21]:
                        final_mutation.append("conservative")
                    else:
                        final_mutation.append("radical")
                final_mut_list.append(final_mutation)

#snp nones

for mutation in point_mutation_list:
    final_mutation = []
    final_mutation += mutation
    count = 0
    for gene in gene_records:
        if int(mutation[4]) > gene[2] and int(mutation[4]) <= gene[3]:
            count+=1
    if count == 0:
        surrounding = mtDNA_seq[int(
            mutation[4]) - 2] + mtDNA_seq[int(mutation[4])]
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append(surrounding)
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mut_list.append(final_mutation)




final_indel_mut_list=[]


#indel Nones

for mutation in indel_mutation_list:
    final_mutation = []
    final_mutation += mutation
    count = 0
    for gene in gene_records:
        if int(mutation[4]) > gene[2] and int(mutation[4]) <= gene[3]:
            count+=1
    if count == 0:
        surrounding = mtDNA_seq[int(
            mutation[4]) - 2] + mtDNA_seq[int(mutation[4])]
        final_mutation.append("None")
        final_mutation.append("None")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append("NA")
        final_mutation.append(len(final_mutation[9])-len(final_mutation[5]))
        final_indel_mut_list.append(final_mutation)


for gene in gene_records:

    #trNA, rRNA, D-loop
    if gene[1] != "CDS":
        gene_length = gene[3] - gene[2]
        for mutation in indel_mutation_list:
            final_mutation = []
            final_mutation += mutation
            mutation_loc = mutation[4]
            if int(mutation_loc) > gene[2] and int(mutation_loc) <= gene[3]:
                ref = int(mutation_loc) - gene[2]
                gene_seq = gene[-1]
                final_mutation.append(str(gene[0]))
                final_mutation.append(str(gene[0]))
                final_mutation.append(
                    str(float(final_mutation[6]) / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(float(final_mutation[6]) / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(1 / region_length_dictionary[gene[1]]))
                final_mutation.append(
                    str(1 / region_length_dictionary[gene[1]]))
                final_mutation.append(len(final_mutation[9])-len(final_mutation[5]))
                final_indel_mut_list.append(final_mutation)

    
# Forward CDS
    else:
        gene_length = gene[3] - gene[2]
        for mutation in indel_mutation_list:
            final_mutation = []
            final_mutation += mutation
            mutation_loc = mutation[4]
            if int(mutation_loc) > gene[2] and int(mutation_loc) <= gene[3]:
                ref = int(mutation_loc) - gene[2]
                # adding polyA tail: needed for some stop codons
                gene_seq = gene[-1] + "AA"
                surrounding = str(gene_seq[ref - 2]) + str(gene_seq[ref])
                b = (ref // 3)
                # Gets the correct codon by using the reference sequence and
                # translates
                right_codon = (gene_seq[(b * 3):(b * 3) + 3])
                codon_list = [right_codon[0], right_codon[1], right_codon[2]]
                # Gets the new(mutated) codon
                codon_list[ref % 3] = mutation[5]
                new_codon = ''
                new_codon = new_codon.join(codon_list)
                original_aa = translate(right_codon.upper())
                mutated_aa = translate(new_codon.upper())
                final_mutation.append(str(gene[0]))
                final_mutation.append(str(gene[1]))
                final_mutation.append(
                    str(float(final_mutation[6]) / gene_length))
                final_mutation.append(str(float(final_mutation[6]) / CDS_len))
                final_mutation.append(str(1 / gene_length))
                final_mutation.append(str(1 / CDS_len))
                final_mutation.append(len(final_mutation[9])-len(final_mutation[5]))
                final_indel_mut_list.append(final_mutation)


for mutation in final_mut_list:
    mutation.append(complement(mutation[13])[1]+complement(mutation[13])[0])

# snp header
mut_header = [
    "old_animal",
    "Animal",
    "tissue",
    "sample_name",
    "ref_num",
    "ref",
    "mut_freq",
    "ref_allele_count",
    "mut_allele_count",
    "mutated_base",
    "som_germ",
    "gene",
    "coding_non",
    "surrounding",
    "mut_freq_div_gene_len",
    "mut_freq_div_len",
    "count_div_gene_len",
    "count_div_len",
    "amino_start",
    "amino_start_group",
    "amino_mut",
    "amino_mut_group",
    "mut_type",
    "conserve_non"]

# Save file
save_file = open("PolG_python_dataset_9_15_20.csv", 'w', newline='')
wr = csv.writer(save_file, quoting=csv.QUOTE_ALL)
wr.writerow(mut_header)
for line in final_mut_list:
    wr.writerow(line)

# indel header
mut_header = [
    "old_animal",
    "Animal",
    "tissue",
    "sample_name",
    "ref_num",
    "ref",
    "mut_freq",
    "ref_allele_count",
    "mut_allele_count",
    "mutated_base",
    "som_germ",
    "i_d",
    "frame",
    "gene",
    "coding_non",
    "mut_freq_div_gene_len",
    "mut_freq_div_len",
    "count_div_gene_len",
    "count_div_len",
    "Change"]

save_file_indels = open(
    "PolG_python_dataset_indels_9_15_20.csv",
    'w',
    newline='')
wr = csv.writer(save_file_indels, quoting=csv.QUOTE_ALL)
wr.writerow(mut_header)
for line in final_indel_mut_list:
    wr.writerow(line)


print("finished!")

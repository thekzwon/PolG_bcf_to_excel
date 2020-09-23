# PolG_bcf_to_excel
This is a python script that takes bcf query files with the fields:

 bcftools query -f '%POS\t%REF\t[%AF]\t[%AD]\t%ALT\n
 
 and outputs 2 excel files, one for point mutations and one for indels.
 
 There should be one query file per sequencing sample with the name:[animal_name]_[Tissue]query.
 
 The column names for the point mutation excel are as follows:
    "old_animal": old name for animal
    "Animal": new animal name, a-z.
    "tissue": Brain or Liver (can be added to)
    "sample_name": [animal_name]+[Tissue]
    "ref_num": number in reference genome
    "ref": reference base
    "mut_freq": mutation frequeny as given by the bcf file
    "ref_allele_count": reference allele count as given by the bcf file
    "mut_allele_count": mutation allele count as given by the bcf file
    "mutated_base": mutated base
    "som_germ":"somatic" or "germline"
    "gene": gene/rRNA/tRNA that the mutation was in, else given NA
    "coding_non": returns CDS/rRNA/tRNA/NA
    "surrounding": gives the bases that surround the mutation
    "mut_freq_div_gene_len": mutation frequency divided by the gene length
    "mut_freq_div_len": mutation frequency divided by the length of the CDS/tRNA/rRNA/
    "count_div_gene_len": mutation count divided by the gene length
    "count_div_len": mutation count divided by the length of the CDS/tRNA/rRNA/
    "amino_start": for CDS regions, gives the starting amino acid
    "amino_start_group": for CDS regions, gives the starting amino acid group, as provided in the supplement
    "amino_mut": for CDS regions, gives the mutated amino acid
    "amino_mut_group":for CDS regions, gives the mutated amino acid group, as provided in the supplement
    "mut_type":silent/missense/nonsense
    "conserve_non":conservative or radical based on the amino acid groups
    
    
The column names for the indel mutation excel are as follows:
    "old_animal":  old name for animal
    "Animal": new animal name, a-z.
    "tissue": Brain or Liver (can be added to)
    "sample_name": [animal_name]+[Tissue]
    "ref_num": number in reference genome
    "ref": reference base
    "mut_freq": mutation frequeny as given by the bcf file
    "ref_allele_count": reference allele count as given by the bcf file
    "mut_allele_count": mutation allele count as given by the bcf file
    "mutated_base": mutated base
    "som_germ": "somatic" or "germline"
    "i_d": insertion or deletion
    "frame":frameshift or not
    "gene": gene/rRNA/tRNA that the mutation was in, else given NA
    "coding_non": returns CDS/rRNA/tRNA/NA
    "mut_freq_div_gene_len": mutation count divided by the gene length
    "mut_freq_div_len": mutation frequency divided by the length of the CDS/tRNA/rRNA/
    "count_div_gene_len": mutation count divided by the length of the CDS/tRNA/rRNA/
    "count_div_len": mutation count divided by the length of the CDS/tRNA/rRNA/
    "Change": number of bases the indel added or deleted
    
    
    
    
    #If there are questions or concerns, please contact Kendra at zwonitz2@utexas.edu
    
    
    
    
    
    
    
    

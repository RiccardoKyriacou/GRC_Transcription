from Bio import SeqIO
from glob import glob
import pandas as pd
import argparse
import os

'''
This script
1) Genrates a table of expressed genes only where expressed genes = genes with a TPM > intergenic_trheshold present in 2 or more samples
2) Parses gene IDs for expressed GRC genes (these are the genes we want to investigate)
3) Using gene ID's parses out the gene sequences from B_cor GRC genome
4) Outputs three FASTA files for GRC genes expressed in soma, germ, and both
5) Parses out all core genes from annotated B_cop genome
6) Produces combined fasta file contining all core chromosome genes
'''
samples_info = {
    # Embryo 0-4h
    "ME1": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME2": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "ME3": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE1": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE2": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    "FE3": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "0-4h"},
    # Embryo 0-8h
    "ML1": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "ML2": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "ML3": {"species": "B_coprophila", "sex": "male", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL1": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL2": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    "FL3": {"species": "B_coprophila", "sex": "female", "tissue": "germ_embryo", "stage": "4-8h"},
    # Late larva/early pupa 
    "Fgerm1": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fgerm2": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fgerm3": {"species": "B_coprophila", "sex": "female", "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Fbody1": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Fbody2": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Fbody3": {"species": "B_coprophila", "sex": "female", "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mgerm1": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mgerm2": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mgerm3": {"species": "B_coprophila", "sex": "male",   "tissue": "germ", "stage": "late-larva-early-pupa"},
    "Mbody1": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mbody2": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    "Mbody3": {"species": "B_coprophila", "sex": "male",   "tissue": "soma_carcass", "stage": "late-larva-early-pupa"},
    # Adult samples for B_coprophila
    "B1":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B2":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B3":  {"species": "B_coprophila", "sex": "male", "tissue": "germ_testes",  "stage": "adult"},
    "B7":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B8":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B9":  {"species": "B_coprophila", "sex": "male", "tissue": "soma_carcass", "stage": "adult"},
    "B10": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B11": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B12": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B16": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B17": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B18": {"species": "B_coprophila", "sex": "gynogenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B19": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B20": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B21": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "germ_ovaries", "stage": "adult"},
    "B25": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B26": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
    "B27": {"species": "B_coprophila", "sex": "androgenic_female", "tissue": "soma_carcass", "stage": "adult"},
}

# 1) Create an output table of expressed genes
# Expressed genes = genes with a TPM > intergenic_trheshold present in 2 or more samples
def get_expression(overlap_file_dir, TPM_threshold):
    expression_files = glob(os.path.join(overlap_file_dir, "gene_overlap*"))
    
    # List to make BLAST queries later on 
    expressed_GRC_gene_lst = [] 

    with open("GRC_gene_expression.tsv", "w") as outf:
        outf.write("gene_id\tscaffold\ttissue_category\tdevelopment_stage\tsex\tgerm_samples\tsoma_samples\tmean_germ_tpm\tmean_soma_tpm\n")
        # Open the output file to write all entries  
        for file in expression_files:
            development_stage = os.path.basename(file).split('_')[2]
            # Open overlap file for given developemnt stage
            with open(file) as f:
                for line in f:
                    columns = line.strip().split("\t")
                    gene_id = columns[0]
                    scaffold = columns[1]
                    sex = columns[3]
                    # Count the number of samples that exceed the threshold
                    count_expresse_germ = 0
                    count_expresse_soma = 0 
                    # Gte TPM
                    sum_germ_TPM = 0
                    sum_soma_TPM = 0 
                    for sample in columns[4:]:
                        if ":" in sample:
                            sample_ID, TPM = sample.split(":")
                            TPM = float(TPM)
                            tissue_exp = samples_info.get(sample_ID, {}).get('tissue', '')

                            if TPM > TPM_threshold and "germ" in tissue_exp:
                                count_expresse_germ += 1
                                sum_germ_TPM += TPM
                            if TPM > TPM_threshold and "soma" in tissue_exp:
                                count_expresse_soma += 1
                                sum_soma_TPM += TPM

                    # Get mean TPMs for samples
                    if count_expresse_germ > 1:
                        mean_germ_TPM = sum_germ_TPM / count_expresse_germ
                    else: mean_germ_TPM = sum_germ_TPM

                    if count_expresse_soma > 1:
                        mean_soma_TPM = sum_soma_TPM / count_expresse_soma
                    else: mean_soma_TPM = sum_soma_TPM

                    # Assign tissue category based on expression rules
                    if count_expresse_germ > 1 and count_expresse_soma == 0:
                        tissue_category = "germ*"
                    elif count_expresse_germ > 1 and count_expresse_soma < 2:
                        tissue_category = "germ"
                    elif count_expresse_soma > 1 and count_expresse_germ < 2:
                        tissue_category = "soma"
                    elif count_expresse_germ > 1 and count_expresse_soma > 1:
                        tissue_category = "both"
                    else:
                        continue  # Skip genes not expressed in at least one category

                    # Write to file
                    outf.write(f"{gene_id}\t{scaffold}\t{tissue_category}\t{development_stage}\t{sex}\t{count_expresse_germ}\t{count_expresse_soma}\t{mean_germ_TPM}\t{mean_soma_TPM}\n")

                    expressed_GRC_gene_lst.append(gene_id)

    return expressed_GRC_gene_lst

# Re-order table 
def sort_expression_table(input_file, output_file):
    # Defin order for development_stage
    stage_order = ["0-4h", "4-8h", "late-larva-early-pupa", "adult"]
    # Load the table into a pandas DataFrame
    df = pd.read_csv(input_file, sep="\t")

    # Ensure `development_stage` is treated as a categorical variable with the defined order
    df["development_stage"] = pd.Categorical(df["development_stage"], categories=stage_order, ordered=True)

    # Sort by gene_id first, then by development_stage
    df_sorted = df.sort_values(by=["gene_id", "development_stage"])

    # Save the sorted table
    df_sorted.to_csv(output_file, sep="\t", index=False)

# 2) Get BLAST quieries for GRC genes and all core genes   
# Function to extract and stitch exons for specified GRC genes
def parse_gff(genome_fasta, gff_file, expressed_GRC_gene_lst):
    # Dictionary to store stitched exon sequences by gene
    GRC_gene_seqs = {}

    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Parse GFF file and extract gene sequences
    with open(gff_file, 'r') as gff:
        for line in gff:
            # Skip any headers/comments
            if line.startswith("#"):
                continue
            
            columns = line.strip().split('\t')
            scaffold = columns[0] 
            feature_type = columns[2]  # Feature type (gene, exon, etc.)
            start = int(columns[3])  # Start position (1-based indexing)
            end = int(columns[4])  # End position
            strand = columns[6]  # Strand (+ or -)

            # Only consider exons
            if feature_type == "gene":
                # Extract the parent gene ID
                gene_ID = columns[8].split("=")[1].split(";")[0]  # Parse ID

                if gene_ID in expressed_GRC_gene_lst:
                    # Get the sequence from the genome
                    gene_sequence = genome[scaffold][start-1:end]  # Adjust to 0-based indexing

                    # Reverse complement if on the negative strand
                    if strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()

                    # Append exonic sequence to the corresponding gene
                    GRC_gene_seqs[gene_ID] = gene_sequence

    return GRC_gene_seqs

# Function to write transcripts to FASTA format
def combined_fasta_output(output_file, transcript_sequences):
    with open(f"{output_file}.fasta", 'w') as out_fasta:
        for gene_id, seq in transcript_sequences.items():
            # Write each transcript in FASTA format
            out_fasta.write(f">{gene_id}\n{str(seq.seq)}\n")

# 2) Functions to generate FASTA file of all core-genome genes
def parse_core_genes_from_gff(genome_fasta, gff_file):
    # Dictionary to store stitched exon sequences by gene
    geneID_sequences = {}

    # Load the genome sequences into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # Parse GFF file and extract exonic sequences
    with open(gff_file, 'r') as gff:
        for line in gff:
            # Skip any headers/comments
            if line.startswith("#"):
                continue
            
            columns = line.strip().split('\t')
            scaffold = columns[0] 
            feature_type = columns[2]  # Feature type (gene, exon, etc.)
            start = int(columns[3])  # Start position (1-based indexing)
            end = int(columns[4])  # End position
            strand = columns[6]  # Strand (+ or -)

            # Only consider full genes
            if feature_type == "gene":
                # Extract the parent gene ID
                gene_ID = columns[8].split("=")[1].split(";")[0]  # Parse ID

                if "s" in gene_ID:
                    # Get the sequence from the genome
                    gene_sequence = genome[scaffold][start-1:end]  # Adjust to 0-based indexing

                    # Reverse complement if on the negative strand
                    if strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()

                    geneID_sequences[gene_ID] = gene_sequence
                        
    return geneID_sequences

def main():
    parser = argparse.ArgumentParser(description="get expressed genes for BLAST")
    parser.add_argument("-g", "--genome", type=str, help="Path to genome to parse", required=True)
    parser.add_argument("-a", "--annotation", type=str, help="Path to annotation", required=True)
    parser.add_argument("-o", "--overlap", type=str, help="Path to dir containing gene_overlap*.tsv files", required=True)
    parser.add_argument("-t", "--threshold", type=float, help="TPM threshold generated from 02_TPM_cutoff", required=False, default = float(0.22))
    args = parser.parse_args()

    # 1) Expression table and GRC gene sequences
    expressed_GRC_genes = get_expression(args.overlap, TPM_threshold=args.threshold)
    sort_expression_table("NEW_GRC_gene_expression.tsv", "NEW_GRC_gene_expression.tsv")
    GRC_gene_seqs = parse_gff(args.genome, args.annotation, expressed_GRC_genes)
    combined_fasta_output("GRC_genes_to_BLAST", GRC_gene_seqs)

    # 2) Core GRC gene fasta
    core_gene_sequences = parse_core_genes_from_gff(args.genome, args.annotation)
    combined_fasta_output("core_genes", core_gene_sequences)


if __name__ == "__main__":
    main()
 
# python 01_get_expressed_genes_for_BLAST.py -g ../../../Annotations/idLycInge5.1.primary.masked.fa -a ../../../Annotations/ling_core_GRC.gff3 -o ../01_STAR_TPM/ 

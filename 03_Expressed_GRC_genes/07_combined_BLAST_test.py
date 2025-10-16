import csv
import argparse

'''
Script that takes in both BLASTp , interpro scan 
Then outputs for each GRC gene of intrest 

Gene_id Scaffold    Development stages  sex core_hit    BlastP_hit  Interpro_hit    

'''

#TODO add repbase results

def parse_BLASTp_output(blast_file, GRC_gene_lst):
    # Store hits
    blast_hits = {}

    with open(blast_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            gene_id = fields[0].split(".")[0]

            # Check for genes of intrest
            if gene_id not in GRC_gene_lst:
                continue
            # eval
            score = float(fields[10])  
            # hit name
            hit_id = fields[14]
            if gene_id not in blast_hits or score > blast_hits[gene_id][1]:
                blast_hits[gene_id] = (hit_id, score)

    return blast_hits


def parse_interpro_output(infile, GRC_gene_lst):
    interpro_hits = {}

    with open(infile) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            gene_id = fields[0].split(".")[0]
            
            if gene_id not in GRC_gene_lst:
                continue

            # Get description 
            ipr_id = fields[12]
            # Get eval 
            score_str = fields[8]  
            if ipr_id != "-" and gene_id not in interpro_hits:
                try:
                    score = float(score_str)
                    interpro_hits[gene_id] = (ipr_id, score)
                except ValueError:
                    continue
    
    return interpro_hits


def write_output(output_file, blast_hits, interpro_hits, GRC_gene_lst):
    # Write final combined summary
    with open(output_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["GENE_ID", "BLASTp_hit", "BLAST_score", "InterPro_hit", "InterPro_score"])

        for gene_id in GRC_gene_lst:
            blast_hit, blast_score = blast_hits.get(gene_id, ("-", "-"))
            interpro_hit, interpro_score = interpro_hits.get(gene_id, ("-", "-"))
            writer.writerow([gene_id, blast_hit, blast_score, interpro_hit, interpro_score])

def main():
    parser = argparse.ArgumentParser(description="get expressed genes for BLAST")
    parser.add_argument("-b", "--blast", type=str, help="Path to blast", required=True)
    parser.add_argument("-i", "--interpro", type=str, help="Path to interpro", required=True)
    args = parser.parse_args()

    output_file = "Combined_BLAST_Interpro_output.tsv"

    GRC_gene_lst = ["g11713","g13362","g13363","g13594","g13694","g14961","g15174","g16029","g17107","g17119","g19121","g19161","g233","g491","g596","g7610","g7957","g8036","g8751"] 
    
    blast_hits = parse_BLASTp_output(args.blast, GRC_gene_lst)
    interpro_hits = parse_interpro_output(args.interpro, GRC_gene_lst)
    write_output(output_file, blast_hits, interpro_hits, GRC_gene_lst)

if __name__ == "__main__":
    main()
# GO_enrichment
GO_enrichment based on the gene annotation by interproscan.

interproscan_GO_Ath.py program


Input files: 
 input_event_genes: file of event-related gene pairs 
 input_interproscan_tsv: tsv result file of interproscan 
 input_lens: lens file of the species 
 go_terms.txt: functional descriptions corresponding to Arabidopsis thaliana GO numbers


Output files: 
 output_1_haveGO_line: line with GO number in tsv file 
 output_2_gene_GO: gene id, corresponding to GO number 
 output_3_GO_gene: GO number, corresponding to gene id 
 output_4_event_genes: pairs of event-related genes, changed to event_genes 
 output_5_event_genes_GO: GO number corresponding to event_genes 
 output_6_event_GO_gene: GO number corresponding to event_genes and descriptive and statistical information 
 output_7_event_genes_GO_gene: further statistical information including Calculation of P-values (hypergeometric distribution)

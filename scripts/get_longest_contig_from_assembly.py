from Bio import SeqIO

assembly_file = snakemake.input.contigs # get the assembly file from from the snakefile input
output_fasta_file = snakemake.output.fasta # get the output path for the fasta file from the snakemake input

longest_contig = None # a counter a variable to hold the longest contig
longest_length = 0 # initiate a counter to store the longest contig length

for record in SeqIO.parse(assembly_file, "fasta"): # for loop to go through each contig in the assembly file
    current_length = len(record.seq) # get the length of the current contig
    if current_length > longest_length: # check if the contig is longer than the current longest contig
        longest_length = current_length # update the longest lenght with the current longest contig
        longest_contig = record # store that contig as the new ongest contig

if longest_contig is not None: # check if the longest contig variable is not empty
    SeqIO.write(longest_contig, output_fasta_file, "fasta") # write the longest contig as a fasta file output

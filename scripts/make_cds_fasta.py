from Bio import SeqIO

cds_fasta=snakemake.input.fasta # get the genome fasta file from from the snakefile input
output_fasta = snakemake.output.cds # get the output path for the fasta file from the snakemake input

with open(output_fasta, "w") as output:
    for record in SeqIO.parse(cds_fasta, "fasta"): # for loop to look through each record in the cds fasta file
        description = record.description # make a variable to hold the the full description line from the fasta header
        if "[protein_id=" in description: # check if the string "[protein_id=" exists in the description
            start = description.find("[protein_id=") + len("[protein_id=") # find the position where "[protein_id=" starts and move the index to just after it
            end = description.find("]", start) # find the closing bracket "]" that marks the end of the protein id
            protein_id = description[start:end] # get the substring between the start and end positions
        else:
            protein_id = record.id  # if no protein id is  found, use the original record ID from the fasta header

        record.id = protein_id # replace the record's ID with the extracted protein id
        record.description = ""  # clear the description so that only the ID appears in the output fasta header
        SeqIO.write(record, output, "fasta") # write the new fasta file with the correct headers

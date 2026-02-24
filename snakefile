FASTQ_DIR = config.get("fastq_dir", "fastq")

def get_samples():
    import pathlib
    p = pathlib.Path(FASTQ_DIR)
    return sorted(f.stem.replace("_1", "") for f in p.glob("*_1.fastq"))

rule all:
    input:
        "PipelineReport.txt",
        "indexes/hcmv_kallisto.idx",
        "results/sleuth_significant.tsv",
        expand("results/kallisto/{sample}", sample=get_samples()),
        expand("counts/{sample}.counts.txt", sample=get_samples()),
        expand("spades/{sample}/contigs.fasta", sample=get_samples()),
        expand("blast/{sample}.blast.tsv", sample=get_samples())

rule download_hcmv_reference:
    output:
        cds="ref/GCF_000845245.1_cds.fna",
        genome="ref/GCF_000845245.1_genome.fna"
    shell:
        """
        mkdir -p ref
        datasets download genome accession GCF_000845245.1 --include genome,cds --filename ref/ncbi_dataset.zip
        unzip -o ref/ncbi_dataset.zip -d ref/
        cp ref/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna {output.cds}
        cp ref/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna {output.genome}
        """

rule get_cds_count:
    input:
        cds="ref/GCF_000845245.1_cds.fna"
    output:
        cds_counts="results/hcmv_cds_counts.txt",
        report="results/report_cds_count.txt"
    shell:
        """
        mkdir -p results
        count=$(grep -c "^>" {input.cds})
        echo "$count" > {output.cds_counts}
        echo "The HCMV genome (GCF_000845245.1) has $count CDS." > {output.report}
        """

rule make_cds_fasta:
    input:
        fasta="ref/GCF_000845245.1_cds.fna"
    output:
        cds="results/hcm_cds.fasta"
    script:
        "scripts/make_cds_fasta.py"

rule kallisto_index:
    input:
        "results/hcm_cds.fasta"
    output:
        "indexes/hcmv_kallisto.idx"
    shell:
        """
        mkdir -p indexes
        kallisto index -i {output} {input}
        """

rule run_kallisto_quant:
    input:
        index="indexes/hcmv_kallisto.idx",
        r1=FASTQ_DIR + "/{sample}_1.fastq",
        r2=FASTQ_DIR + "/{sample}_2.fastq"
    output:
        directory("results/kallisto/{sample}")
    shell:
        """
        kallisto quant -i {input.index} -b 30 -o {output} {input.r1} {input.r2}
        """

rule run_sleuth:
    input:
        expand("results/kallisto/{sample}", sample=get_samples())
    output:
        "results/sleuth_significant.tsv"
    script:
        "scripts/sleuth.R"

rule bowtie2_index:
    input:
        "ref/GCF_000845245.1_genome.fna"
    output:
        directory("indexes/bowtie2_hcmv")
    shell:
        """
        mkdir -p {output}
        bowtie2-build {input} {output}/bowtie2_index
        """

rule bowtie2_map:
    input:
        r1=FASTQ_DIR + "/{sample}_1.fastq",
        r2=FASTQ_DIR + "/{sample}_2.fastq",
        index="indexes/bowtie2_hcmv"
    output:
        r1="results/bowtie2/{sample}.mapped_1.fastq",
        r2="results/bowtie2/{sample}.mapped_2.fastq"
    shell:
        """
        mkdir -p results/bowtie2
        bowtie2 -x {input.index}/bowtie2_index -1 {input.r1} -2 {input.r2} --no-unal | samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n
        """

rule count_reads:
    input:
        r1=FASTQ_DIR + "/{sample}_1.fastq",
        mapped_r1="results/bowtie2/{sample}.mapped_1.fastq"
    output:
        counts="counts/{sample}.counts.txt",
        report="results/report_reads_{sample}.txt"
    shell:
        """
        mkdir -p counts
        original=$(($(wc -l < {input.r1})/4))
        mapped=$(($(wc -l < {input.mapped_r1})/4))
        echo "original_reads: $original" > {output.counts}
        echo "mapped_reads: $mapped" >> {output.counts}
        echo "Sample {wildcards.sample} had $original read pairs before and $mapped read pairs after Bowtie2 filtering" > {output.report}
        """

rule spades_assembly:
    input:
        r1="results/bowtie2/{sample}.mapped_1.fastq",
        r2="results/bowtie2/{sample}.mapped_2.fastq"
    output:
        "spades/{sample}/contigs.fasta"
    shell:
        """
        mkdir -p spades/{wildcards.sample}
        spades.py -1 {input.r1} -2 {input.r2} --only-assembler -k 127 -o spades/{wildcards.sample}
        """

rule get_longest_contig:
    input:
        contigs="spades/{sample}/contigs.fasta"
    output:
        fasta="blast/{sample}.longest_contig.fasta"
    script:
        "scripts/get_longest_contig_from_assembly.py"

rule make_blast_db:
    output:
        "db/betaherpesvirinae.nsq"
    shell:
        """
        mkdir -p db
        datasets download virus genome taxon betaherpesvirinae --refseq --include genome --filename db/db.zip
        unzip -o db/db.zip -d db
        makeblastdb -in db/ncbi_dataset/data/*.fna -dbtype nucl -out db/betaherpesvirinae
        """

rule blast_longest_contig:
    input:
        fasta="blast/{sample}.longest_contig.fasta",
        db="db/betaherpesvirinae.nsq"
    output:
        "blast/{sample}.blast.tsv"
    shell:
        """
        mkdir -p blast
        blastn -query {input.fasta} -db db/betaherpesvirinae -out {output} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 5 -max_hsps 1
        """

rule assemble_pipeline_report:
    input:
        "results/report_cds_count.txt",
        expand("results/report_reads_{sample}.txt", sample=get_samples()),
        "results/sleuth_significant.tsv",
        expand("blast/{sample}.blast.tsv", sample=get_samples())
    output:
        "PipelineReport.txt"
    shell:
        """
        cat results/report_cds_count.txt > {output}
        echo -e "\ntarget_id\ttest_stat\tpval\tqval" >> {output}
        tail -n +2 results/sleuth_significant.tsv >> {output}
        echo "" >> {output}
        cat results/report_reads_*.txt >> {output}
        for f in blast/*.blast.tsv; do
            sample=$(basename $f .blast.tsv)
            echo -e "\n$sample:" >> {output}
            echo -e "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle" >> {output}
            head -n 5 $f >> {output}
        done
        """

rule cleanup:
    shell:
        """
        rm -rf ncbi_dataset ncbi_dataset.zip
        rm -rf ref
        rm -rf indexes
        rm -rf results
        rm -rf counts
        rm -rf spades
        rm -rf blast
        rm -rf db
        rm -f PipelineReport.txt
        """

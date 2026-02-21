SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

rule all:
    input:
        "PipelineReport.txt",
        "indexes/hcmv_kallisto.idx",
        "results/sleuth_significant.tsv",
        expand("results/kallisto/{sample}", sample=SAMPLES),
        expand("counts/{sample}.counts.txt", sample=SAMPLES),
        expand("spades/{sample}/contigs.fasta", sample=SAMPLES),
        expand("blast/{sample}.blast.tsv", sample=SAMPLES)

rule download_hcmv_reference:
    output:
        cds="ref/GCF_000845245.1_cds.fna",
        genome="ref/GCF_000845245.1_genome.fna"
    shell:
        """
        mkdir -p ref/GCF_000845245.1
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
        r1="fastq/{sample}_1.fastq",
        r2="fastq/{sample}_2.fastq"
    output:
        directory("results/kallisto/{sample}")
    shell:
        """
        kallisto quant -i {input.index} -b 30 -o {output} {input.r1} {input.r2}
        """

rule run_sleuth:
    input:
        expand("results/kallisto/{sample}", sample=SAMPLES)
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
        r1="fastq/{sample}_1.fastq",
        r2="fastq/{sample}_2.fastq",
        index="indexes/bowtie2_hcmv"
    output:
        "results/bowtie2/{sample}.mapped.bam"
    shell:
        """
        mkdir -p results/bowtie2
        bowtie2 -x {input.index}/bowtie2_index -1 {input.r1} -2 {input.r2} --no-unal | samtools view -bS - > {output}
        """

rule count_reads:
    input:
        r1="fastq/{sample}_1.fastq",
        bam="results/bowtie2/{sample}.mapped.bam"
    output:
        counts="counts/{sample}.counts.txt",
        report="results/report_reads_{sample}.txt"
    shell:
        """
        mkdir -p counts
        original=$(($(wc -l < {input.r1})/4))
        mapped=$(($(samtools view -c {input.bam}) / 2))
        echo "original_reads: $original" > {output.counts}
        echo "mapped_reads: $mapped" >> {output.counts}
        echo "Sample {wildcards.sample} had $original read pairs before and $mapped read pairs after Bowtie2 filtering" > {output.report}
        """

rule spades_assembly:
    input:
        "results/bowtie2/{sample}.mapped.bam"
    output:
        "spades/{sample}/contigs.fasta"
    shell:
        """
        mkdir -p spades/{wildcards.sample}
        samtools fastq -1 spades/{wildcards.sample}/r1.fq -2 spades/{wildcards.sample}/r2.fq -0 /dev/null -s /dev/null {input}
        spades.py -1 spades/{wildcards.sample}/r1.fq -2 spades/{wildcards.sample}/r2.fq --only-assembler -k 127 -o spades/{wildcards.sample}
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
        blastn -query {input.fasta} -db db/betaherpesvirinae -out {output} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 5
        """

rule assemble_pipeline_report:
    input:
        "results/report_cds_count.txt",
        expand("results/report_reads_{sample}.txt", sample=SAMPLES),
        "results/sleuth_significant.tsv",
        expand("blast/{sample}.blast.tsv", sample=SAMPLES)
    output:
        "PipelineReport.txt"
    shell:
        """
        cat results/report_cds_count.txt > {output}
        echo "" >> {output}
        cat results/report_reads_*.txt >> {output}
        echo -e "\ntarget_id\ttest_stat\tpval\tqval" >> {output}
        tail -n +2 results/sleuth_significant.tsv >> {output}
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
        rm -rf indexes/bowtie2_hcmv
        rm -rf results/bowtie2/*
        rm -rf counts/*
        rm -rf results/kallisto/*
        rm -rf spades/*
        rm -rf blast/*
        rm -f results/hcm_cds.fasta
        rm -f results/hcmv_cds_counts.txt
        rm -f indexes/hcmv_kallisto.idx
        rm -rf db/*
        """

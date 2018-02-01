SAMPLES = ["SRR1635435", "SRR1635436", "SRR1635459", "SRR1635460"]

rule all:
	input:
		expand("fastqc/{sample}_1_fastqc.html", sample=SAMPLES),
		expand("bow/{sample}_1_unique.bam", sample=SAMPLES)

rule sra_to_fastq:
	input:
		"data/{sample}"
	output:
		"fastq/{sample}_1.fastq.gz"
		# temp("fastq/{sample}_2.fastq.gz")
	log:
		"logs/fastq-dump/{sample}.log"
	params:
		options="--outdir fastq --gzip --skip-technical --readids --origfmt --dumpbase --split-files --clip"
	shell:
		"fastq-dump {params.options} {input} &> {log}"

rule fastq_to_fastqc:
	input:
		"fastq/{sample}_1.fastq.gz"
	output:
		"fastqc/{sample}_1_fastqc.html"
	shell:
		"fastqc -o fastqc {input}"

rule trimgalore:
	input:
		"fastq/{sample}_1.fastq.gz"
	output:
		"trim_galore/{sample}_1_trimmed.fq.gz"
	shell:
		"trim_galore --length 16 --output_dir trim_galore {input}"

rule bowtie2_index:
	input:
		"genome/human.fa.gz"
	output:
		"bow/human.*.bt2"
	shell:
		"bowtie2-build {input} bow/human"

rule bowtie2_align:
	input:
		"trim_galore/{sample}_1_trimmed.fq.gz"
	output:
		samfile="bow/{sample}_1_unique.sam",
		unique="bow/{sample}_1_unique.bam",
		output="bow/{sample}_1_output.bam",
		header="bow/{sample}_1_header.bam"
	threads: 2
	shell:
		"""
		bowtie2 -x human/human -S {output.samfile} -p {threads} -U {input}
		samtools view -bS -t human/human.fa.fai {output.samfile} > {output.output}
		samtools view -H {output.output} > {output.header}
		samtools view -F 4 {output.output} | grep -v "XS:" | cat {output.header} - | samtools view -b - > {output.unique}
		rm {output.header}
		"""

SAMPLES = ["SRR1635435", "SRR1635436", "SRR1635459", "SRR1635460"]

rule all:
	input:
		expand("fastqc/{sample}_1_fastqc.html", sample=SAMPLES),
		expand("trim_galore/{sample}_1_trimmed.fq.gz", sample=SAMPLES)

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

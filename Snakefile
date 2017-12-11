SAMPLES = ["SRR1635435", "SRR1635436", "SRR1635459", "SRR1635460"]

rule all:
	input:
		expand("fastq/{sample}_1.fastq.gz", sample=SAMPLES),
		# expand("fastq/{sample}_2.fastq.gz", sample=SAMPLES)

rule sra_to_fastq:
	input:
		"data/{sample}"
	output:
		temp("fastq/{sample}_1.fastq.gz"),
		# temp("fastq/{sample}_2.fastq.gz")
	params:
		fd="--outdir fastq --gzip --skip-technical --readids --origfmt --dumpbase --split-files --clip"
	shell:
		"fastq-dump {params.fd} {input}"

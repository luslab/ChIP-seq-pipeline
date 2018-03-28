import os
import re
from itertools import chain

# Import config file & parameters
configfile: "config.yaml"

# Import paths from config file
PATH = "{}".format(config['global']['snakemake'])
RESULTS_PATH = "{}".format(config['global']['results'])
GENOME_PATH = "{}/{}".format(config['genome']['path'], config['genome']['name'])
ANNOTATION_PATH = "{}/{}".format(config['annotation']['path'], config['annotation']['name'])
BT2_INDEX_PATH = "{}/{}".format(config['alignment']['path'], config['alignment']['bowtie2'])

rule all:
	input:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/nrf/{sample}_1_NRF.txt", sample=config["samples"]),
		RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak",
		expand(RESULTS_PATH+"/fastq/{sample}_1.fastq.gz", sample=config["samples"])

rule map:
	input:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"])

# Fastq files
rule sra_to_fastq:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
	log:
		RESULTS_PATH+"/logs/fastq-dump/{sample}.log"
	params:
		default="--outdir {}/fastq {}".format(RESULTS_PATH, config['parameters']['fastqc'])
	run:
		if len(re.findall("([.fq|.fastq]+[.gz]*)$", input[0])) == 0:
			# if os.path.isfile(input[0]):
			shell("fastq-dump {params.default} {input} &> {log}")
			# else:
			# 	srafile = os.path.basename(input[0])
			# 	shell("fastq-dump {params.default} {srafile} &> {log}")
		else:
			shell("ln -s {} {} &> {}".format(os.path.abspath(input[0]), os.path.abspath(output[0]), log[0]))

# Quality control
rule fastq_to_fastqc:
	input:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
	output:
		RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html",
		temp(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.zip")
	params:
		default="-o {}/fastqc {}".format(RESULTS_PATH, config['parameters']['fastq_dump'])
	shell:
		"fastqc {params.default} {input}"

# Adapter removal
rule trimgalore:
	input:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
	output:
		trimmed=temp(RESULTS_PATH+"/trimmed/{sample}_1_trimmed.fq.gz"),
		report=RESULTS_PATH+"/trimmed/{sample}_report.txt"
	log:
		RESULTS_PATH+"/logs/trim_galore/{sample}.log"
	params:
		default="--output_dir {}/trimmed {}".format(RESULTS_PATH, config['parameters']['trim_galore'])
	run:
		shell("trim_galore {params.default} {input} &> {log}")
		report_name="/".join(( os.path.split(output['trimmed'])[0], os.path.basename(input[0])+"_trimming_report.txt"))
		shell("mv {report_name} {output.report}")

# Genome index
rule bowtie2_index:
	input:
		GENOME_PATH
	output:
		protected(expand(BT2_INDEX_PATH+".{count}.bt2", count=['1', '2', '3', '4'])),
		protected(expand(BT2_INDEX_PATH+".rev.{count}.bt2", count=['1', '2']))
	run:
		shell("bowtie2-build {input} "+BT2_INDEX_PATH)

# Alignment
rule bowtie2_align:
	input:
		RESULTS_PATH+"/trimmed/{sample}_1_trimmed.fq.gz"
	output:
		samfile=temp(RESULTS_PATH+"/mapped/{sample}_1.sam"),
		sorted=temp(RESULTS_PATH+"/mapped/{sample}_1.bam"),
		unique=RESULTS_PATH+"/mapped/{sample}_1_unique.bam",
		index=RESULTS_PATH+"/mapped/{sample}_1_unique.bam.bai"
	threads: 8
	log:
		RESULTS_PATH+"/logs/bowtie2/{sample}.log"		
	params:
		default="".format(config['parameters']['bowtie2'])
	run:
		shell("bowtie2 -x {BT2_INDEX_PATH} -S {output.samfile} -p {threads} {params.default} -U {input} 2> {log}")
		shell("samtools view -b {output.samfile} | samtools sort -O BAM -o {output.sorted} -")
		shell("samtools view -F 4 -h {output.sorted} | grep -E '@|NM:' | grep -v 'XS:' | samtools view -b > {output.unique}")
		shell("samtools index {output.unique}")

# Non-redundant fraction
rule NRF:
	input:
		RESULTS_PATH+"/mapped/{sample}_1_unique.bam"
	output:
		RESULTS_PATH+"/nrf/{sample}_1_NRF.txt"
	script:
		PATH+"/scipts/NRF.py"

# Peak calling
rule peak_calling:
	input:
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"])
	output:
		RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak",
		RESULTS_PATH+"/peaks/macs2_peaks.xls",
		RESULTS_PATH+"/peaks/macs2_summits.bed",
		RESULTS_PATH+"/peaks/macs2_model.r"
	log:
		RESULTS_PATH+"/logs/macs2/peak_calling.log"
	params:
		default="--outdir {}/peaks {}".format(RESULTS_PATH, config['parameters']['macs2'])
	run:
		control = list(filter(lambda x: 'Control' in x, input))
		target = list(filter(lambda x: 'Control' not in x, input))
		shell('macs2 callpeak -t {target} -c {control} -f BAM {params.default} -n macs2 > {log}')

# Individual dag
rule dag:
	output:
		RESULTS_PATH+"/dag.pdf"
	params:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/nrf/{sample}_1_NRF.txt", sample=config["samples"]),
	run:
		shell("snakemake --dag {} | dot -Tpdf > {}".format(" ".join([ item[0] for item in params]), *output))

# Complete dag
rule dag_complete:
	output:
		RESULTS_PATH+"/dag_complete.pdf"
	params:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/nrf/{sample}_1_NRF.txt", sample=config["samples"]),
		[RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak"]
	run:
		shell("snakemake --dag {} | dot -Tpdf > {}".format(" ".join(list(chain.from_iterable(params))), *output))

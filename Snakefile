from itertools import chain

# Import config file & parameters
configfile: "config.yaml"

# Import paths from config file
# PATH = config['snakemake_path']
PATH = "{}".format(config['global']['snakemake'])
# RESULTS_PATH = config['results_path']
RESULTS_PATH = "{}".format(config['global']['results'])
# GENOME_PATH = config['genome_path']
GENOME_PATH = "{}/{}".format(config['genome']['path'], config['genome']['name'])
# ANNOTATION_PATH = config['annotation_path']
ANNOTATION_PATH = "{}/{}".format(config['annotation']['path'], config['annotation']['name'])
# STAR_INDEX_PATH = config['star_index_path']
BT2_INDEX_PATH = "{}/{}".format(config['alignment']['path'], config['alignment']['bowtie2'])

rule all:
	input:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/nrf/{sample}_1_NRF.txt", sample=config["samples"]),
		RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak"

rule sra_to_fastq:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
		# temp("fastq/{sample}_2.fastq.gz")
	log:
		RESULTS_PATH+"/logs/fastq-dump/{sample}.log"
	params:
		options="--outdir fastq --gzip --skip-technical --readids --origfmt --dumpbase --split-files --clip"
	shell:
		"fastq-dump {params.options} {input} &> {log}"

rule fastq_to_fastqc:
	input:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
	output:
		RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html",
		temp(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.zip")
	shell:
		"fastqc -o fastqc {input}"

rule trimgalore:
	input:
		RESULTS_PATH+"/fastq/{sample}_1.fastq.gz"
	output:
		RESULTS_PATH+"/trimmed/{sample}_1_trimmed.fq.gz"
	log:
		RESULTS_PATH+"/logs/trim_galore/{sample}.log"
	shell:
		"""
		trim_galore --length 16 --output_dir trim_galore {input}
		mv {input}_trimming_report.txt {log}
		"""

rule bowtie2_index:
	input:
		GENOME_PATH
	output:
		expand(BT2_INDEX_PATH+".{count}.bt2", count=['1', '2', '3', '4']),
		expand(BT2_INDEX_PATH+".rev.{count}.bt2", count=['1', '2'])
	run:
		shell("bowtie2-build {input} "+BT2_INDEX_PATH)

rule bowtie2_align:
	input:
		RESULTS_PATH+"/trimmed/{sample}_1_trimmed.fq.gz"
	output:
		samfile=temp(RESULTS_PATH+"/mapped/{sample}_1_unique.sam"),
		unique=RESULTS_PATH+"/mapped/{sample}_1_unique.bam",
		output=temp(RESULTS_PATH+"/mapped/{sample}_1_output.bam"),
		header=temp(RESULTS_PATH+"/mapped/{sample}_1_header.bam")
	threads: 8
	shell:
		"""
		bowtie2 -x human/human -S {output.samfile} -p {threads} -U {input}
		samtools view -bS -t human/human.fa.fai {output.samfile} > {output.output}
		samtools view -H {output.output} > {output.header}
		samtools view -F 4 {output.output} | grep -v "XS:" | cat {output.header} - | samtools view -b - > {output.unique}
		"""
rule NRF:
	input:
		RESULTS_PATH+"/mapped/{sample}_1_unique.bam"
	output:
		RESULTS_PATH+"/nrf/{sample}_1_NRF.txt"
	script:
		PATH+"/NRF.py"

rule peak_calling:
	input:
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"])
	output:
		RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak"
	params:
		outdir=RESULTS_PATH
	run:
		control = list(filter(lambda x: 'Control' in x, input))
		target = list(filter(lambda x: 'Control' not in x, input))
		shell('macs2 callpeak -t {target} -c {control} -f BAM --outdir {params.outdir} -n macs2')

rule dag:
	output:
		RESULTS_PATH+"/dag.pdf"
	params:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_NRF.txt", sample=config["samples"]),
	run:
		shell("snakemake --dag {} | dot -Tpdf > {}".format(" ".join([ item[0] for item in params]), *output))

rule dag_complete:
	output:
		RESULTS_PATH+"/dag_complete.pdf"
	params:
		expand(RESULTS_PATH+"/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_unique.bam", sample=config["samples"]),
		expand(RESULTS_PATH+"/mapped/{sample}_1_NRF.txt", sample=config["samples"]),
		[RESULTS_PATH+"/peaks/macs2_peaks.narrowPeak"]
	run:
		shell("snakemake --dag {} | dot -Tpdf > {}".format(" ".join(list(chain.from_iterable(params))), *output))

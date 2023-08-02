
# Run ____________________ snakemake -s quality_control.smk --use-conda 

###############################################################################
### CONFIGURATION ###
#####################

configfile: "config.yaml"

###############################################################################
### FUNCTIONS ###
#################

###############################################################################
### WILDCARDS ###
#################
FASTQ = glob_wildcards("resources/reads/{fastq}.fastq.gz")

###############################################################################
### RESOURCES ###
#################

# OS = config["os"]                  # Operating system
CPUS = config["resources"]["cpus"] # Threads (maximum)
# REFERENCE = config["consensus"]["reference"] # Reference genome

###############################################################################
### ENVIRONMENTS ###
####################

NANOPLOT = config["conda"]["osx"]["nanoplot"]         # NanoPlot conda env


###############################################################################
### RULES ###
#############

rule all:
    input: expand("results/00_Quality_Control/{fastq}/{fastq}_NanoStats.txt", fastq = FASTQ.fastq)
       
rule nanoplot:
    # Aim: Generate quality control plots for nanopore reads
    # Use: NanoPlot --fastq <fastq> --threads <threads> --outdir <outdir> --prefix <prefix>
    message: "NanoPlot - Quality control plots for nanopore reads"
    conda:
        NANOPLOT
    input:
        fastq = "resources/reads/{fastq}.fastq.gz"
        # bam = "results/02_Mapping/{reference}/{fastq}_mark-dup.bam"
    output:
        stats = "results/00_Quality_Control/{fastq}/{fastq}_NanoStats.txt"
    log:
        "results/10_Reports/tools-log/nanoplot/{fastq}_nanoplot.log"
    shell:
        "NanoPlot "
        "--fastq {input.fastq} "
        # "--bam {input.bam} "
        "--threads {CPUS} "
        "--outdir results/00_Quality_Control/{wildcards.fastq}/ "
        "--prefix {wildcards.fastq}_ "
        "&> {log}"
    

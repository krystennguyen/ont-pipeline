###I###RAA###D######U###2###3###3#######T###R###A###N###S###V###I###H###M###I####
###                                                                         ###
###    /\  ______      ___ ____ _  _ __   ____ __   ____     ______  /\     ###
###    ||  \ \ \ \    / __| ___| \/ )__\ (  _ (  ) (_  _)   / / / /  ||     ###
###    ||   > > > >  ( (_-.)__) \  /(__)\ )   /)(__ _)(_   < < < <   ||     ###
###    ||  /_/_/_/    \___(____) \(__)(__|_)\_|____|____)   \_\_\_\  ||     ###
###    \/                                                            \/     ###
###                                                                         ###
###I###R###D######U###2###3###3#######T###R###A###N###S###V###I###H###M###I####
# Name ___________________ quality_control.smk
# Version ________________ v.2023.06
# Author _________________ Nicolas Fernandez
# Affiliation ____________ IRD_U233_TransVIHMI
# Aim ____________________ Snakefile with quality control rules
# Date ___________________ 2021.09.28
# Latest modifications ___ 2023.06.21
# Run ____________________ snakemake -s quality_control.smk --use-conda 

###############################################################################
### CONFIGURATION ###
#####################

configfile: "configuration/config.yaml"

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

OS = config["os"]                  # Operating system
CPUS = config["resources"]["cpus"] # Threads (maximum)

###############################################################################
### ENVIRONMENTS ###
####################

MULTIQC = config["conda"][OS]["multiqc"]           # Multi-QC conda env
FASTQ_SCREEN = config["conda"][OS]["fastq_screen"] # Fastq-Screen conda env
FASTQC= config["conda"][OS]["fastqc"]              # FastQC conda env

###############################################################################
### PARAMETERS ###
##################

SUBSET = config["fastq_screen"]["subset"]     # Fastq-Screen --subset
FQC_CONFIG = config["fastq_screen"]["config"] # Fastq-Screen --conf
#MQC_CONFIG = config["multiqc"]["config"]      # MultiQC --conf
#TAG = config["multiqc"]["tag"]                # MultiQC --tag

###############################################################################
### RULES ###
#############

rule all:
    input:
        multiqc = "results/00_Quality_Control/multiqc/",
        fastq_screen = expand("results/00_Quality_Control/fastq-screen/{fastq}/",
                             fastq = FASTQ),
        fastqc = expand("results/00_Quality_Control/fastqc/{fastq}/",
                        fastq = FASTQ)

###############################################################################
rule multiqc_reports_aggregation:
    # Aim: aggregates bioinformatics analyses results into a single report
    # Use: multiqc [OPTIONS] --output [MULTIQC/] [FASTQC/] [MULTIQC/]
    message:
        "MultiQC reports aggregating"
    conda:
        MULTIQC
    #params:
        #config = MQC_CONFIG,
        #tag = TAG
    input:
        fastqc = expand("results/00_Quality_Control/fastqc/{fastq}/",
                        fastq = FASTQ),
        fastq_screen = expand("results/00_Quality_Control/fastq-screen/{fastq}/",
                             fastq = FASTQ)
    output:
        multiqc = directory("results/00_Quality_Control/multiqc/")
    log:
        "results/10_Reports/tools-log/multiqc.log"
    shell:
        "multiqc "                  # Multiqc, searches in given directories for analysis & compiles a HTML report
        "--quiet "                   # -q: Only show log warning
        "--no-ansi "                 # Disable coloured log
        #"--config {params.config} "  # Specific config file to load
        #"--tag {params.tag} "        # Use only modules which tagged with this keyword
        #"--pdf "                     # Creates PDF report with 'simple' template (Requires Pandoc to be installed)
        #"--export "                  # Export plots as static images in addition to the report
        "--outdir {output.multiqc} " # -o: Create report in the specified output directory
        "{input.fastqc} "            # Input FastQC files
        "{input.fastq_screen} "      # Input Fastq-Screen
        "&> {log}"                   # Log redirection

###############################################################################
rule fastqscreen_contamination_checking:
    # Aim: screen if the composition of the library matches with  what you expect
    # Use fastq_screen [OPTIONS] --outdir [DIR/] [FASTQ.GZ]
    message:
        "Fastq-Screen contamination check for [ {wildcards.fastq} ] reads"
    conda:
        FASTQ_SCREEN
    resources:
        cpus = CPUS
    params:
        config = FQC_CONFIG,
        subset = SUBSET
    input:
        fastq = "resources/reads/{fastq}.fastq.gz"
    output:
        fastq_screen = directory("results/00_Quality_Control/fastq-screen/{fastq}/")
    log:
        "results/10_Reports/tools-log/fastq-screen/{fastq}.log"
    shell:
        "fastq_screen "                  # FastqScreen, what did you expect ?
        "-q "                             # --quiet: Only show log warning
        "--threads {resources.cpus} "     # --threads: Specifies across how many threads bowtie will be allowed to run
        "--aligner 'bwa' "                # -a: choose aligner 'bowtie', 'bowtie2', 'bwa'
        "--conf {params.config} "         # path to configuration file
        "--subset {params.subset} "       # Don't use the whole sequence file, but create a subset of specified size
        "--outdir {output.fastq_screen} " # Output directory
        "{input.fastq} "                  # Input file.fastq
        "&> {log}"                        # Log redirection

###############################################################################
rule fastqc_quality_control:
    # Aim: reads sequence files and produces a quality control report
    # Use: fastqc [OPTIONS] --output [DIR/] [FASTQ.GZ]
    message:
        "FastQC quality control for [ {wildcards.fastq} ] reads"
    conda:
        FASTQC
    resources:
        cpus = CPUS
    input:
        fastq = "resources/reads/{fastq}.fastq.gz"
    output:
        fastqc = directory("results/00_Quality_Control/fastqc/{fastq}/")
    log:
        "results/10_Reports/tools-log/fastqc/{fastq}.log"
    shell:
        "mkdir -p {output.fastqc} "    # (*) this directory must exist as the program will not create it
        "2> /dev/null && "             # in silence and then... 
        "fastqc "                    # FastQC, a high throughput sequence QC analysis tool
        "--quiet "                    # -q: Supress all progress messages on stdout and only report errors
        "--threads {resources.cpus} " # -t: Specifies files number which can be processed simultaneously
        "--outdir {output.fastqc} "   # -o: Create all output files in the specified output directory (*)
        "{input.fastq} "              # Input file.fastq
        "&> {log}"                    # Log redirection

###############################################################################

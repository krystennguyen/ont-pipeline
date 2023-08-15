###############################################################################
### CONFIGURATION ###
#####################

configfile: "config.yaml"

###############################################################################
### FUNCTIONS ###
#################

def get_memory_per_thread(wildcards):
    memory_per_thread = int(RAM) // int(CPUS)
    return memory_per_thread

###############################################################################
### WILDCARDS ###
#################
SAMPLE = glob_wildcards("resources/reads/{sample}.fastq.gz")


###############################################################################
### RESOURCES ###
#################

OS = config["conda"]["osx"]              # Operating system
CPUS = config["resources"]["cpus"]       # Threads (maximum)
RAM = config["resources"]["ram"]         # Memory (RAM) in Gb (maximum)
MEM_GB = get_memory_per_thread           # Memory per thread in GB (maximum)
TMP_DIR = config["resources"]["tmp_dir"] # Temporary directory

###############################################################################
### ENVIRONMENTS ###
####################

# PYCHOPPER = config["conda"]["osx"]["pychopper"] # Pychopper conda environment
MINIMAP2 = config["conda"]["osx"]["minimap2"]
SAMTOOLS = config["conda"]["osx"]["samtools"]
BEDTOOLS = config["conda"]["osx"]["bedtools"]
GAWK = config["conda"]["osx"]["gawk"]
# MEDAKA = config["conda"]["osx"]["medaka"]
NANOPOLISH = config["conda"]["osx"]["nanopolish"]


###############################################################################
### PARAMETERS ###
##################
PRIMER = ("resources/primers/fasta/PCS111_primers.fas")  # Path to primer
REF_PATH = config["consensus"]["path"]          # Path to genomes references
REFERENCE = config["consensus"]["reference"]    # Genome reference sequence, in fasta format
MIN_COV = config["consensus"]["min_cov"]        # Minimum coverage, mask lower regions with 'N' 
MODEL = config["consensus"]["model"]            # Model for consensus sequence calling


###############################################################################
### RULES ###
#############
rule all:
    input:
        masked_ref = expand("results/04_Variants/{reference}/{sample}_{min_cov}X_masked-ref.fasta", sample = SAMPLE.sample, reference = REFERENCE, min_cov = MIN_COV),
        cov_stats = expand("results/03_Coverage/{reference}/{sample}_{min_cov}X_coverage-stats.tsv", sample = SAMPLE.sample, reference = REFERENCE, min_cov = MIN_COV),
        consensus = expand("results/05_Consensus/{reference}/{sample}_{min_cov}X_consensus.fasta", sample = SAMPLE.sample, reference = REFERENCE, min_cov = MIN_COV)
###############################################################################
rule nanopolish_variant:
    # Aim: 
    # Use: nanopolish variants [OPTIONS] --reads [READS.fasta] --bam [ALIGNMENTS.bam] --genome [REFERENCE.fasta]
    message:
        "Nanopolish variant calling for [[ {wildcards.sample} ]] sample (for {wildcards.reference}, @{wildcards.min_cov}X)"
    conda: 
        NANOPOLISH
    params:
        ref = expand("{ref_path}{reference}.fasta", ref_path = REF_PATH, reference = REFERENCE)
    input:
        reads = "results/04_Variants/{reference}/{sample}_{min_cov}X_masked-ref.fasta",
        bam = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
    output:
        consensus = "results/05_Consensus/{reference}/{sample}_{min_cov}X_consensus.fasta"
    log:
        "results/10_Reports/tools-log/medaka/{reference}/{sample}_{min_cov}X_variant.log"
    shell:
        "nanopolish variant "
        "--consensus "
        "--reads {input.reads} "
        "--bam {input.bam} "
        "--genome {params.ref} "
        "&> {log}"  
    


# rule medaka_variant:
#     # Aim: variant calling
#     # Use: medaka variant [REFERENCE.fasta] [CONSENSUS.hdf] [OUTPUT.vcf]
#     message:
#         "Medaka variant calling for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
#     conda:
#         MEDAKA
#     params:
#         ref = expand("{ref_path}{reference}.fasta", ref_path = REF_PATH, reference = REFERENCE)
#     input:
#         consensus = "results/05_Consensus/{reference}/{sample}_consensus.hdf"
#     output:
#         variant = "results/04_Variants/{reference}/{sample}_variant.vcf"
#     log:
#         "results/10_Reports/tools-log/medaka/{reference}/{sample}_variant.log"
#     shell:
#         "medaka variant " # Medaka variant, tools for variant calling
#         "--check_output "        # Check output file
#         "{params.ref} " # Reference fasta file
#         "{input.consensus} " # Consensus hdf input
#         "{output.variant} " # Variant output
#         "&> {log}"                    # Log redirection
# ###############################################################################
# rule medaka_consensus:
#     # Aim: consensus sequence
#     # Use: medaka consensus -i [MAPPED.bam] -d [REFERENCE.fasta] -o [CONSENSUS.fasta]
#     message:
#         "Medaka consensus sequence for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
#     conda:
#         MEDAKA
#     params:
#         ref = expand("{ref_path}{reference}.fasta", ref_path = REF_PATH, reference = REFERENCE),
#         model = MODEL
#     input:
#         mapped = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
#     output:
#         consensus = "results/05_Consensus/{reference}/{sample}_consensus.hdf"
#     log:
#         "results/10_Reports/tools-log/medaka/{reference}/{sample}_consensus.log"
#     shell:
#         "medaka consensus " # Medaka consensus, tools for consensus sequence calling
#         "--debug "
#         # "--threads 2 "                  # Number of threads
#         # "-d {params.ref} "   # Reference fasta file
#         # "--check_output "        # Check output file
#         "--model {params.model} " # Model for consensus sequence calling
#         "--chunk_len 800 "     # Chunk length
#         "--chunk_ovlp 400 "        # Chunk overlap
#         "{input.mapped} " # Mapped bam input
#         "{output.consensus} " # Consensus output
#         "&> {log}"                    # Log redirection


###############################################################################
rule bedtools_masking:
    # Aim: masking low coverage regions
    # Use: bedtools maskfasta [OPTIONS] -fi [REFERENCE.fasta] -bed [RANGE.bed] -fo [MASKEDREF.fasta]
    message:
        "BedTools masking low coverage regions for [[ {wildcards.sample} ]] sample (for {wildcards.reference}, @{wildcards.min_cov}X)"
    conda:
        BEDTOOLS
    params:
        path = REF_PATH
    input:
        low_cov_mask = "results/03_Coverage/{reference}/bed/{sample}_{min_cov}X_low-cov-mask.bed"
    output:
        masked_ref = "results/04_Variants/{reference}/{sample}_{min_cov}X_masked-ref.fasta"
    log:
        "results/10_Reports/tools-log/bedtools/{reference}/{sample}_{min_cov}X_masking.log"
    shell:
        "bedtools maskfasta "                       # Bedtools maskfasta, mask a fasta file based on feature coordinates
        "-fi {params.path}{wildcards.reference}.fasta " # Input FASTA file 
        "-bed {input.low_cov_mask} "                 # BED/GFF/VCF file of ranges to mask in -fi
        "-fo {output.masked_ref} "                   # Output masked FASTA file
        "&> {log}"                                   # Log redirection 

###############################################################################
rule bedtools_merged_mask:
    # Aim: merging overlaps
    # Use: bedtools merge [OPTIONS] -i [FILTERED.bed] -g [GENOME.fasta] 
    message:
        "BedTools merging overlaps for [[ {wildcards.sample} ]] sample (for {wildcards.reference}, @{wildcards.min_cov}X)"
    conda:
        BEDTOOLS
    input:
        min_cov_filt = "results/03_Coverage/{reference}/bed/{sample}_{min_cov}X_min-cov-filt.bed"
    output:
        low_cov_mask = temp("results/03_Coverage/{reference}/bed/{sample}_{min_cov}X_low-cov-mask.bed")
    log:
        "results/10_Reports/tools-log/bedtools/{reference}/{sample}_{min_cov}X_merging.log"
    shell:
        "bedtools merge "          # Bedtools merge, merges overlapping BED/GFF/VCF entries into a single interval
        "-i {input.min_cov_filt} "  # -i: BED/GFF/VCF input to merge 
        "1> {output.low_cov_mask} " # merged output
        "2> {log}"                  # Log redirection

###############################################################################
rule awk_min_covfilt:
    # Aim: minimum coverage filtration
    # Use: awk '$4 < [MIN_COV]' [BEDGRAPH.bed] 1> [FILTERED.bed]
    message:
        "Awk minimum coverage filtration for [[ {wildcards.sample} ]] sample (for {wildcards.reference}, @{wildcards.min_cov}X)"
    conda:
        GAWK
    input:
        genome_cov = "results/03_Coverage/{reference}/bed/{sample}_genome-cov.bed"
    output:
        min_cov_filt = temp("results/03_Coverage/{reference}/bed/{sample}_{min_cov}X_min-cov-filt.bed")
    log:
        "results/10_Reports/tools-log/awk/{reference}/{sample}_{min_cov}X_min-cov-filt.log"
    shell:
        "awk "                      # Awk, a program that you can use to select particular records in a file and perform operations upon them
        "'$4 < {wildcards.min_cov}' " # Minimum coverage for masking regions in consensus sequence
        "{input.genome_cov} "         # BedGraph coverage input
        "1> {output.min_cov_filt} "   # Minimum coverage filtered bed output
        "2> {log} "                   # Log redirection

###############################################################################
rule awk_coverage_statistics:
    # Aim: computing genomme coverage stats
    # Use: awk {FORMULA} END {{print [RESULTS.tsv] [BEDGRAPH.bed]
    message:
        "Awk compute genome coverage statistics BED [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
    conda:
        GAWK
    input:
        #sickle = "results/10_Reports/tools-log/sickle-trim/{sample}.log",
        samtools = "results/10_Reports/tools-log/samtools/{reference}/{sample}_mark-dup.log",
        flagstat = "results/03_Coverage/{reference}/flagstat/{sample}_flagstat.json",
        histogram = "results/03_Coverage/{reference}/histogram/{sample}_coverage-histogram.txt",
        genome_cov = "results/03_Coverage/{reference}/bed/{sample}_genome-cov.bed"
    output:
        cov_stats = "results/03_Coverage/{reference}/{sample}_{min_cov}X_coverage-stats.tsv"
    log:
        "results/10_Reports/tools-log/awk/{reference}/{sample}_{min_cov}X_coverage-stats.log"
    shell:
        # """ rawReads=$(grep -o -E  """                                  # Get raw reads 
        # """ 'Total read pairs processed:.+' {input.cutadapt}  """       #
        # """ | sed -E 's/Total read pairs processed:\ +//'  """          #
        # """ | sed 's/,//g') ; """                                       #
        #
        # """ cutadaptPF=$(grep -o -E """                                 # Get cutadapt Passing Filtes reads
        # """ 'Pairs written \(passing filters\):.+' {input.cutadapt} """ #
        # """ | sed -E 's/Pairs written \(passing filters\):\ +//' """    #
        # """ | sed 's/,//g') ; """                                       #
        #
        # """ sicklePF=$(grep -o -E """                                   # Get sickle Passing Filtes reads
        # """ 'FastQ paired records kept:.+' {input.sickle} """           #
        # """ | sed -E 's/FastQ paired records kept:\ +//') ; """         #
        #
        """ totalDuplicate=$(grep -o -E """                             # Get total duplicated reads
        """ 'DUPLICATE TOTAL:.+' {input.samtools} """                   #
        """ | sed -E 's/DUPLICATE TOTAL:\ +//') ; """                   #
        #
        """ estimatedLibrarySize=$(grep -o -E """                       # Get estimated library size
        """ 'ESTIMATED_LIBRARY_SIZE:.+' {input.samtools} """            #
        """ | sed -E 's/ESTIMATED_LIBRARY_SIZE:\ +//') ; """            #
        #
        """ samtoolsPF=$(grep -o -E """                                 # Get samtool Passing Filter reads
        """ 'WRITTEN: .+' {input.samtools} """                          #
        """ | sed -E 's/WRITTEN:\ +//') ; """                           #
        #
        """ mappedReads=$(grep -o -E -m 1 """                           # Get mapped reads
        """ '"mapped": .+' {input.flagstat} """                         #
        """ | sed -E 's/"mapped":\ +//' """                             #
        """ | sed 's/,//g') ; """                                       #
        #
        """ mappedPercentReads=$(grep -o -E -m 1 """                    # Get mapped precent reads
        """ '"mapped %": .+' {input.flagstat} """                       #
        """ | sed -E 's/"mapped %":\ +//' """                           #
        """ | sed 's/,//g') ; """                                       #
        #
        """ covPercentAt1X=$(grep -o -E """                             # Get coverage percent @1X
        """ 'Percent covered:.+' {input.histogram} """                  #
        """ | sed -E 's/Percent covered:\ +//') ; """                   #
        #
        """ awk """                                                   # Awk, a program to select particular records in a file and perform operations upon them
        # """ -v rawReads="${{rawReads}}" """                             # Define external variable
        # """ -v cutadaptPF="${{cutadaptPF}}" """                         # Define external variable
        # """ -v sicklePF="${{sicklePF}}" """                             # Define external variable
        """ -v totalDuplicate="${{totalDuplicate}}" """                 # Define external variable
        """ -v estimatedLibrarySize="${{estimatedLibrarySize}}" """     # Define external variable
        """ -v samtoolsPF="${{samtoolsPF}}" """                         # Define external variable
        """ -v mappedReads="${{mappedReads}}" """                       # Define external variable
        """ -v mappedPercentReads="${{mappedPercentReads}}" """         # Define external variable
        """ -v covPercentAt1X="${{covPercentAt1X}}" """                 # Define external variable
        """ '$4 >= {wildcards.min_cov} {{supMin_Cov+=$3-$2}} ; """      # Genome size (>= min_cov @X)
        """ {{genomeSize+=$3-$2}} ; """                                 # Genome size (total)
        """ {{totalBases+=($3-$2)*$4}} ; """                            # Total bases @1X
        """ {{totalBasesSq+=(($3-$2)*$4)**2}} """                       # Total bases square @1X
        """ END """                                                    # END
        """ {{print """                                                # Print
        """ "sample_id", "\t", """                                      # header: Sample ID
        """ "raw_paired_reads", "\t", """                               # header: Raw paired reads
        # """ "cutadapt_pairs_PF", "\t", """                              # header: Cutadapt Passing Filters
        # """ "sickle_reads_PF", "\t", """                                # header: Sickle-trim Passing Filters
        """ "duplicated_reads", "\t", """                               # header:
        # """ "duplicated_percent_%","\t", """                            # header:
        """ "estimated_library_size*", "\t", """                        # header:
        """ "samtools_pairs_PF", "\t", """                              # header:
        # """ "mapped_with", "\t",  """                                   # header: aligner
        """ "mapped_on", "\t",  """                                     # header: reference
        """ "mapped_reads", "\t", """                                   # header:
        """ "mapped_percent_%", "\t", """                               # header:
        """ "mean_depth", "\t", """                                     # header: mean depth
        """ "standard_deviation", "\t", """                             # header: standard deviation
        """ "cov_percent_%_@1X", "\t", """                              # header: coverage percentage @1X
        """ "cov_percent_%" "\t", """                                   # header: coverage percentage
        """ "@_min_cov" """                                             # header: @_[min_cov]_X
        """ ORS """                                                      # \n newline
        """ "{wildcards.sample}", "\t", """                             # value: Sample ID
        """ rawReads, "\t", """                                         # value: Raw sequences
        # """ cutadaptPF, "\t", """                                       # value: Cutadapt Passing Filter
        # """ sicklePF, "\t", """                                         # value: Sickle Passing Filter
        """ totalDuplicate, "\t", """                                   # value:
        # """ int(((totalDuplicate)/(rawReads*2))*100), "%", "\t", """    # value: (divided by 2 to estimated pairs)
        """ estimatedLibrarySize, "\t", """                             # value:
        """ samtoolsPF, "\t", """                                       # value:
        # """ "{wildcards.aligner}", "\t",  """                           # value: aligner
        """ "{wildcards.reference}", "\t",  """                         # value: reference
        """ mappedReads, "\t", """                                      # value:
        """ mappedPercentReads, "%", "\t", """                          # value:
        """ int(totalBases/genomeSize), "\t", """                       # value: mean depth
        """ int(sqrt((totalBasesSq/genomeSize)-(totalBases/genomeSize)**2)), "\t", """ # Standard deviation value
        """ covPercentAt1X, "\t", """                                   # value
        """ supMin_Cov/genomeSize*100, "%", "\t", """                   # Coverage percent (@ min_cov X) value
        """ "@{wildcards.min_cov}X" """                                 # @ min_cov X value
        """ }}' """                                                     # Close print
        """ {input.genome_cov} """                                      # BedGraph coverage input
        """ 1> {output.cov_stats} """                                   # Mean depth output
        """ 2> {log}"""                                                 # Log redirection
        
###############################################################################
rule bedtools_genome_coverage:
    # Aim: computing genome coverage sequencing
    # Use: bedtools genomecov [OPTIONS] -ibam [MARK-DUP.bam] 1> [BEDGRAPH.bed]
    message:
        "BedTools computing genome coverage for [[ {wildcards.sample} ]] sample against reference genome sequence (for {wildcards.reference})"
    conda:
        BEDTOOLS
    input:
        mark_dup = "results/02_Mapping/{reference}/{sample}_mark-dup.bam",
        index = "results/02_Mapping/{reference}/{sample}_mark-dup.bam.bai"
    output:
        genome_cov = "results/03_Coverage/{reference}/bed/{sample}_genome-cov.bed"
    log:
        "results/10_Reports/tools-log/bedtools/{reference}/{sample}_genome-cov.log"
    shell:
        "bedtools genomecov "    # Bedtools genomecov, compute the coverage of a feature file among a genome
        "-bga "                   # Report depth in BedGraph format, regions with zero coverage are also reported
        "-ibam {input.mark_dup} " # The input file is in BAM format, must be sorted by position
        "1> {output.genome_cov} " # BedGraph output
        "2> {log} "               # Log redirection
###############################################################################
rule samtools_coverage_histogram:
    # Aim: alignment depth and percent coverage histogram
    # Use: samtools coverage --histogram [INPUT.bam]
    message:
        "SamTools calcul alignment depth and percent coverage @1X from BAM file for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
       #bins = BINS,
       #depth = DEPTH
    input:
        mark_dup = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
    output:
        histogram = "results/03_Coverage/{reference}/histogram/{sample}_coverage-histogram.txt"
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_coverage-histogram.log"
    shell:
        "samtools coverage "          # Samtools coverage, tools for alignments in the SAM format with command to alignment depth and percent coverage
        "--histogram "                 # -m: show histogram instead of tabular output
        "--verbosity 4 "               # Set level of verbosity [INT] (default: 3)
        "--n-bins 149 "                # -w: number of bins in histogram (default: terminal width - 40) (todo: {params.bins}) 
        "--depth 0 "                   # -d maximum allowed coverage depth [INT] (default: 1000000 ; 0 removing any depth limit) (todo: {params.depth}) 
        "--output {output.histogram} " # write output to FILE (default: stdout)
        "{input.mark_dup} ; "          # Mark_dup bam input
        "echo >> {output.histogram} ; " # Newline
        "samtools coverage "    # Samtools coverage, tools for alignments in the SAM format with command to alignment depth and percent coverage
        "--verbosity 4 "         # Set level of verbosity [INT] (default: 3)
        "{input.mark_dup} "      # Mark_dup bam input
        ">> {output.histogram} " # write output to FILE (default: stdout)
        "2> {log}"               # Log redirection

###############################################################################
rule samtools_flagstat_ext:
    # Aim: simple stats
    # Use: samtools flagstat -@ [THREADS] [INPUT.bam]
    message:
        "SamTools calcul simple stats from BAM file for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        mark_dup = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
    output:
        flagstat = "results/03_Coverage/{reference}/flagstat/{sample}_flagstat.{ext}"
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_flagstat-{ext}.log"
    shell:
        "samtools flagstat "           # Samtools flagstat, tools for alignments in the SAM format with command to simple stat
        "--threads {resources.cpus} "   # -@: Number of additional threads to use (default: 1)
        "--verbosity 4 "                # Set level of verbosity [INT] (default: 3)
        "--output-fmt {wildcards.ext} " # -O Specify output format (none, tsv and json)
        "{input.mark_dup} "             # Mark_dup bam input
        "1> {output.flagstat} "        # Mark_dup index output
        "2> {log}"                      # Log redirection

###############################################################################
rule samtools_index_markdup:
    # Aim: indexing marked as duplicate BAM file
    # Use: samtools index -@ [THREADS] -b [MARK-DUP.bam] [INDEX.bai]
    message:
        "SamTools indexing marked as duplicate BAM file for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        mark_dup = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
    output:
        index = "results/02_Mapping/{reference}/{sample}_mark-dup.bam.bai"
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_mark-dup-index.log"
    shell:
        "samtools index "     # Samtools index, tools for alignments in the SAM format with command to index alignment
        "-@ {resources.cpus} " #--threads: Number of additional threads to use (default: 1)(NB, --threads form dose'nt work)
        "-b "                  # -b: Generate BAI-format index for BAM files (default)
        "{input.mark_dup} "    # Mark_dup bam input
        "{output.index} "      # Mark_dup index output
        "&> {log}"             # Log redirection

###############################################################################
rule samtools_markdup:
    # Aim: marking duplicate alignments
    # Use: samtools markdup -@ [THREADS] -r -s -O BAM [SORTED.bam] [MARK-DUP.bam] 
    message:
        "SamTools marking duplicate alignments for [[ {wildcards.sample} ]] sample (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        sorted = "results/02_Mapping/{reference}/{sample}_sorted.bam"
    output:
        mark_dup = "results/02_Mapping/{reference}/{sample}_mark-dup.bam"
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_mark-dup.log"
    shell:
        "samtools markdup "          # Samtools markdup, tools for alignments in the SAM format with command mark duplicates
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-r "                         # -r: Remove duplicate reads
        "-s "                         # -s: Report stats
        "--output-fmt BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.sorted} "             # Sorted bam input
        "{output.mark_dup} "          # Mark_dup bam output
        "&> {log}"                    # Log redirection 

###############################################################################
rule samtools_sorting:
    # Aim: sorting
    # Use: samtools sort -@ [THREADS] -m [MEM_GB] -T [TMP_DIR] -O BAM -o [SORTED.bam] [FIX-MATE.bam] 
    message:
        "SamTools sorting [[ {wildcards.sample} ]] sample reads (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
       mem_gb = MEM_GB,
       tmp_dir = TMP_DIR
    input:
        fix_mate = "results/02_Mapping/{reference}/{sample}_fix-mate.bam"
    output:
        sorted = temp("results/02_Mapping/{reference}/{sample}_sorted.bam")
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_sorted.log"
    shell:
        "samtools sort "             # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "     # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-T {resources.tmp_dir} "     # -T: Write temporary files to PREFIX.nnnn.bam
        "-O BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sorted} "         # Sorted bam output
        "{input.fix_mate} "           # Fixmate bam input
        "&> {log}"                    # Log redirection 

###############################################################################
rule samtools_fixmate:
    # Aim: filling in mate coordinates
    # Use: samtools fixmate -@ [THREADS] -m -O BAM [SORT-BY-NAMES.bam] [FIX-MATE.bam] 
    message:
        "SamTools filling in mate coordinates [[ {wildcards.sample} ]] sample reads (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS
    input:
        sort_by_names = "results/02_Mapping/{reference}/{sample}_sort-by-names.bam"
    output:
        fix_mate = temp("results/02_Mapping/{reference}/{sample}_fix-mate.bam")
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_fix-mate.log"
    shell:
        "samtools fixmate "          # Samtools fixmate, tools for alignments in the SAM format with command to fix mate information
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-m "                         # -m: Add mate score tag 
        "-O BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "{input.sort_by_names} "      # Sort_by_names bam input
        "{output.fix_mate} "          # Fix_mate bam output 
        "&> {log}"                    # Log redirection 

###############################################################################
rule samtools_sortbynames:
    # Aim: sorting by names
    # Use: samtools sort -t [THREADS] -m [MEM_GB] -n -O BAM -o [SORT-BY-NAMES.bam] [MAPPED.sam]
    message:
        "SamTools sorting by names [[ {wildcards.sample} ]] sample reads (for {wildcards.reference})"
    conda:
        SAMTOOLS
    resources:
       cpus = CPUS,
       mem_gb = MEM_GB
    input:
        mapped = "results/02_Mapping/{reference}/{sample}_mapped.sam"
    output:
        sort_by_names = temp("results/02_Mapping/{reference}/{sample}_sort-by-names.bam")
    log:
        "results/10_Reports/tools-log/samtools/{reference}/{sample}_sort-by-names.log"
    shell:
        "samtools sort "             # Samtools sort, tools for alignments in the SAM format with command to sort alignment file
        "--threads {resources.cpus} " # -@: Number of additional threads to use (default: 1)
        "-m {resources.mem_gb}G "     # -m: Set maximum memory per thread, suffix K/M/G recognized (default: 768M)
        "-n "                         # -n: Sort by read name (not compatible with samtools index command) 
        "--output-fmt BAM "           # -O: Specify output format: SAM, BAM, CRAM (here, BAM format)
        "-o {output.sort_by_names} "  # -o: Write final output to FILE rather than standard output
        "{input.mapped} "             # Mapped reads input
        "&> {log}"                    # Log redirection 

###############################################################################
rule minimap2:
    # Aim: Align sequences against reference data
    # Use: minimap2 -x map-ont -a [REFERENCE.fasta] [INPUT.fastq.gz] > [OUTPUT.sam]
    message:
        "Minimap2 mapping sequence against reference data"
    conda:
        MINIMAP2
    params:
        ref= expand("{ref_path}{reference}.fasta", ref_path = REF_PATH, reference = REFERENCE)
    input:
        reads = "resources/reads/{sample}.fastq.gz"
    output:
        mapped = ("results/02_Mapping/{reference}/{sample}_mapped.sam")
    log:
        "results/10_Reports/tools-log/minimap2/{reference}/{sample}_mapped.log"
    shell:
        "minimap2 "
        "-x map-ont "       # Preset Nanopore reference mapping
        "-a "               # Output in SAM format
        "{params.ref} "     # Reference fasta file
        "{input.reads} "     # Reads input
        "-o {output.mapped} " # Mapped output
        "&> {log}"                    # Log redirection 


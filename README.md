# Analysis pipeline for nanopore sequence data
Input:
~/resources/reads -- fastq.gz files of basecalled, demultiplexed reads
~/resources/primers -- primer bed/bedpe/fasta files
~/resources/genomes -- reference genome fasta files

Ouput: ~/results
</br>

Use: 
```
bash start_pipeline.sh
```
</br>

Features:
- Quality control with NanoPlot
- Reads mapping with Minimap2
- Genome coverage with Samtools and Bedtools
- Variant calling with Nanopolish (TODO)
- Consensus sequences Nanopolish (TODO)
- Lineage reports (TODO)




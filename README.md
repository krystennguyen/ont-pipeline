# Analysis pipeline for nanopore sequence data
Input:
</br>
~/resources/reads -- fastq.gz files of basecalled, demultiplexed reads </br>
~/resources/signal_data -- raw Nanopore signal data fast5 files (Nanopolish workflow)
~/resources/primers -- primer bed/bedpe/fasta files </br>
~/resources/genomes -- reference genome fasta files </br>


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
- Variant calling (TODO)
- Consensus sequences(TODO)
- Lineage reports (TODO)

Bugs to fix:
- Medaka illegal hardware instructions
- 




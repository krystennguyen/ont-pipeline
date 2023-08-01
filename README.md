# Analysis pipeline for nanopore sequence data
Input: resources/reads  ----- fastq.gz files of demultiplexed reads
</br>
Use:
```
snakemake --cores 8 -s workflow/snakefiles/pipeline.smk --use-conda  
```

Features:
- Reads mapping
- Genome coverage
- Variant calling
- Consensus sequences
- Lineage reports


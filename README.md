# Benchmarking GEM-X and PIPseq V
Here, we benchmark 10x Genomics GEM-X 3' scRNA-seq against Illumina/PIP-seq V
using PBMCs from the same donor. These are cursory analyses in the
Seurat environment ot assess UMI abundance and patterns of gene expression. 

Overall:
- Cell type abundances, annotated by azimuth, are overall consistent.
- 10x Genomics GEM-X has a pretty substantial UMI detection advantage
- The pipelines yield overall concordant gene expression results (r = 0.89) but many genes are selectively 
	- Not entirely clear yet if this is algorithmic or chemistry associated 

## Reproducibility

- [Download large data files from dropbox here](https://www.dropbox.com/scl/fo/ggzq4kkuat11sj4n5p4x9/AEHxXG7ASikewV2KbQ8kzFY?rlkey=scc0doxmqvubof17drs195n84&dl=0)
- Organize folders like so

```
├── benchmark_large_data (dropbox contents)
└── singlecell_genomics_benchmark (this repository)
```

## Upstream processing

### PIPseq v5
```
pipseeker-v3.3.0-linux/pipseeker full --chemistry V --fastq KK-3734_PIPseq_library --star-index-path pipseeker-gex-reference-GRCh38-2022.04 --output-path PIPseq_v5_50 --threads 50
```

### 10x Genomics
```
cr="~/cellranger-8.0.1/cellranger"
tr="~/references/txg/refdata-gex-GRCh38-2024-A"
s="KK-3733_PIPseq_comparison_GEMX_IGO_12437_JK_3"

$cr count --fastqs=fastqs --id Comparison_3p --transcriptome $tr --create-bam false --sample $s
```

## Questions? Comments
Contact: [lareauc@mskcc.org](lareauc@mskcc.org)
<br>
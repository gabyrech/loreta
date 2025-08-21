# loreta — ONT targeted long‑read analysis \


Pipeline for analysis of Oxford Nanopore long reads in targeted assays (e.g. Cas9 targeting, Adaptive Sampling).


## Highlights
- Adapter trimming (Porechop) and quality filtering (fastp)
- Read length histograms (full, middle 80%, bottom 80%)
- NanoPlot QC on filtered reads
- Minimap2 mapping to reference and SOI (sequence of interest)
- Coverage computation (mosdepth) and BigWig generation (deepTools)
- PAF alignments for Ref/SOI
- Extraction of reads mapping to SOI and re‑mapping those reads to the reference
- SV calling with Sniffles2 on SOI and SOI→Ref BAMs
- Event detection from per‑base coverage at user‑supplied thresholds with gene overlap annotation
- Optional mapping + SV calling to an additional reference (e.g. CHM13) if provided
- Summary report mirroring the original `00_summary.txt`


## Quick start
```bash
# run with conda
nextflow run . -profile conda -params-file params.json


# run with Singularity/Apptainer (HPC)
nextflow run . -profile singularity,slurm -params-file params.json


# run with Docker (workstation)
nextflow run . -profile docker -params-file params.json
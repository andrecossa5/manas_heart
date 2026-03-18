# nf-manas-heart

Nextflow pipeline for somatic variant calling from LCM-WGS data (human heart).
Merges BAMs by chunk (heart) or cell type (placenta trophoblasts), then runs Mutect2 on each heart chunk against the merged placenta normal.

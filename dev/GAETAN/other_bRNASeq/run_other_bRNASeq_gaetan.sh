#nextflow run nf-core/rnaseq -r 3.19.0 -name rnaseq_gaetan2 -profile docker -work-dir /home/mna_bioinfo/Bureau/work/MAIN/runrnaseq_gaetan2 -params-file /media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/GAETAN/other_bRNASeq/nf-params_other_bRNASeq.json

nextflow run nf-core/rnaseq -r 3.19.0 resume rnaseq_gaetan2 -profile docker -work-dir /home/mna_bioinfo/Bureau/work/MAIN/runrnaseq_gaetan2 -params-file /media/mna_bioinfo/MNA2_Stockage/CHAZAUD/dev/GAETAN/other_bRNASeq/nf-params_other_bRNASeq.json

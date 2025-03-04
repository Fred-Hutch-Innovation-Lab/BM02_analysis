#!/bin/bash
SLURM_OPTS="--job-name=satcurvepipseeker \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=8 \
             --mem-per-cpu=8G \
             --time=40:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"

samples=(B2-FBv4_F1A B2-FBv4_F1B B2-FBv4_F5A B2-FBv4_F5B)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)
sample_cell_counts=(["B2-FBv4_F1A"]=18953 ["B2-FBv4_F1B"]=16574 ["B2-FBv4_F5A"]=16646 ["B2-FBv4_F5B"]=15151)
total_reads=(["B2-FBv4_F1A"]=449738139  ["B2-FBv4_F1B"]=415600000 ["B2-FBv4_F5A"]=419575000 ["B2-FBv4_F5B"]=380850000)
for sample in ${samples[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/fluent_v4/sat_curves/${sample}"
    mkdir $outdir;
    cd $outdir;
    for target in ${targets[@]}; do
        reads=$((target * ${sample_cell_counts[$sample]}));
        sbatch $SLURM_OPTS \
        --wrap "/fh/fast/_IRC/FHIL/grp/bioinfo_tools/pipseeker/pipseeker-v3.3.0-linux/pipseeker full \
                --fastq /fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/fluent_v4/downsampled/${sample} \
                --star-index-path /shared/biodata/ngs/Reference/PIPseeker/pipseeker-gex-reference-GRCh38-2022.04 \
                --chemistry v4 \
                --annotation /shared/biodata/ngs/Reference/PIPseeker/human-pbmc-references/references/human-pbmc-v4-detailed.csv \
                --output-path ${outdir}/${target} \
                --remove-bam \
                --force-cells ${sample_cell_counts[$sample]} \
                --random-seed 33 \
                --input-reads ${total_reads[$sample]} \
                --downsample-to ${reads}";
    done
done

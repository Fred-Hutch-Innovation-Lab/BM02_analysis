#!/bin/bash
SLURM_OPTS="--job-name=satcurvegemx \
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

samples=(GEMX_B23XF1A GEMX_B23XF1B GEMX_B23XF5A GEMX_B23XF5B)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)
sample_cell_counts=(["GEMX_B23XF1A"]=37147  ["GEMX_B23XF1B"]=38360 ["GEMX_B23XF5A"]=32383 ["GEMX_B23XF5B"]=31155)
total_reads=(["GEMX_B23XF1A"]=945800000  ["GEMX_B23XF1B"]=968075000 ["GEMX_B23XF5A"]=823200000 ["GEMX_B23XF5B"]=786625000)
for sample in ${samples[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/gemx_3p/sat_curves/${sample}"
    # mkdir $outdir;
    for target in ${targets[@]}; do
        outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/gemx_3p/sat_curves/${sample}/${target}"
        # mkdir $outdir;
        cd $outdir;
        reads=$(echo "scale=5; $target * ${sample_cell_counts[$sample]} / ${total_reads[$sample]}" | bc);
        sed -e "s/{{SAMPLE}}/${sample:9:3}/g" -e "s/{{SUBSAMPLE}}/$reads/g" /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/gemx_3p/sat_curves/config.csv > $outdir/config.csv
        sbatch $SLURM_OPTS \
        --wrap "cellranger multi \
            --id=$sample --csv=$outdir/config.csv";
    done
done

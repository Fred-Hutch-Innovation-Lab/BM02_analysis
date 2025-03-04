#!/bin/bash
module load CellRanger/8.0.0

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

samples=(F1A F1B F5A F5B)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)

declare -A sample_cell_counts 
sample_cell_counts=(["F1A"]=16338 ["F1B"]=16403 ["F5A"]=13078 ["F5B"]=13413)
# total_reads=(["F1A"]=411350000 ["F1B"]=412625000 ["F5A"]=328950000 ["F5B"]=336875000)
for sample in ${samples[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/nextgem_3p/sat_curves/${sample}"
    mkdir $outdir;
    for target in ${targets[@]}; do
        outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/nextgem_3p/sat_curves/${sample}/${target}"
        mkdir $outdir;
        cd $outdir;
        # reads=$(echo "scale=5; $target * ${sample_cell_counts[$sample]}" | bc);
        # echo ${sample_cell_counts[${sample}]};
        sed -e "s/{{SAMPLE}}/${sample}/g" -e "s/{{TARGET}}/$target/g" -e "s/{{CELLS}}/${sample_cell_counts[$sample]}/g" /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/nextgem_3p/sat_curves/config_manual_downsample.csv > $outdir/config.csv
        sbatch $SLURM_OPTS \
        --wrap "cellranger multi \
            --id=$sample --csv=$outdir/config.csv";
    done
done

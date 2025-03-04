#!/bin/bash
# sleep 4h
SLURM_OPTS="--job-name=satcurveflex \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=16 \
             --mem-per-cpu=8G \
             --time=24:00:00 
             "

ml CellRanger/8.0.0
samples=(F1A F1B F5A F5B)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000)
# sample_cell_counts=(["F1A"]=16338 ["F1B"]=16403 ["F5A"]=13078 ["F5B"]=13413)
# total_reads=(["F1A"]=411350000 ["F1B"]=412625000 ["F5A"]=328950000 ["F5B"]=336875000)
for target in ${targets[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/flex/sat_curves/${target}"
    mkdir $outdir;
    cd $outdir;
    # reads=$(echo "scale=5; $target * ${sample_cell_counts[$sample]} / ${total_reads[$sample]}" | bc);
        sed -e "s/{{TARGET}}/$target/g" /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/flex/sat_curves/config.csv > $outdir/config.csv
        sbatch $SLURM_OPTS \
        --wrap "cellranger multi --csv=$outdir/config.csv --id=flex_sat_${target}";
done

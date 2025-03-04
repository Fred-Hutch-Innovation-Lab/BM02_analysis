#!/bin/bash
SLURM_OPTS="--job-name=satcurvegemx \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=8 \
             --mem-per-cpu=8G \
             --time=10:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"
module load CellRanger/8.0.0
samples=(F1 F2 F4 F5)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)
declare -A sample_cell_counts 
sample_cell_counts=(["F1"]=18103 ["F2"]=22708 ["F4"]=21211 ["F5"]=19601)
total_reads=(["F1"]=468500000  ["F2"]=578825000 ["F4"]=538300000 ["F5"]=490400000)
for sample in ${samples[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/gemx_5p/sat_curves/${sample}"
    mkdir $outdir;
    for target in ${targets[@]}; do
        outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/gemx_5p/sat_curves/${sample}/${target}"
        mkdir $outdir;
        cd $outdir;
        # reads=$(echo "scale=5; $target * ${sample_cell_counts[$sample]} / ${total_reads[$sample]}" | bc);
        sed -e "s/{{SAMPLE}}/${sample}/g" -e "s/{{TARGET}}/$target/g" -e "s/{{CELLS}}/${sample_cell_counts[$sample]}/g" /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/gemx_5p/sat_curves/config.csv > $outdir/config.csv
        sbatch $SLURM_OPTS \
        --wrap "cellranger multi \
            --id=$sample --csv=$outdir/config.csv";
    done
done

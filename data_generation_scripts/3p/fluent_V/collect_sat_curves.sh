
indir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/fluent_V/sat_curves"
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)

echo "kit,sample,target,umis,genes" > $indir/fluent_V_satcurves.csv
for directory in $indir/B2-FBv5_F*; do
    sample=$(basename $directory)
    for target in ${targets[@]}; do
        results="$directory/$target/metrics/force_*/metrics_summary.csv";
        umi=$(grep "median_transcripts_in_cells" $results | \
        # awk -F'"' '
        # {
        #     for (i = 2; i <= NF; i += 2) {
        #         gsub(",", "", $i)
        #     }
        #     print
        # }' | \
        cut -d ',' -f 2) #| \
        #

        genes=$(grep "median_genes_in_cells" $results | \
        # awk -F'"' '
        # {
        #     for (i = 2; i <= NF; i += 2) {
        #         gsub(",", "", $i)
        #     }
        #     print
        # }' | cut -d ',' -f 6 |
        # awk '{ gsub(/^[ \t]+|[ \t]+$/, "", $0); print }');
        cut -d ',' -f 2) 
        echo "Fluent_V,$sample,$target,$umi,$genes" >> $indir/fluent_V_satcurves.csv
    done
done
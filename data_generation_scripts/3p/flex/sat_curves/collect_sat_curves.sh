indir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/flex/sat_curves"
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000)

echo "kit,sample,target,genes,umis" > $indir/flex_3p_satcurves.csv
for sample in F1A F1B F5A F5B; do
    # sample=$(basename $directory)
    for target in ${targets[@]}; do
        results=$indir/$target/flex_sat_$target/outs/per_sample_outs/$sample/metrics_summary.csv;
        umi=$(grep "Median UMI counts per cell" $results | \
        awk -F'"' '
        {
            for (i = 2; i <= NF; i += 2) {
                gsub(",", "", $i)
            }
            print
        }' | cut -d ',' -f 6 |
        awk '{ gsub(/^[ \t]+|[ \t]+$/, "", $0); print }');

        genes=$(grep "Median genes per cell" $results | \
        awk -F'"' '
        {
            for (i = 2; i <= NF; i += 2) {
                gsub(",", "", $i)
            }
            print
        }' | cut -d ',' -f 6 |
        awk '{ gsub(/^[ \t]+|[ \t]+$/, "", $0); print }');
        echo "Flex,Flex_$sample,$target,$genes,$umi" >> $indir/flex_3p_satcurves.csv
    done
done

indir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/nextgem_3p/sat_curves"
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)

echo "kit,sample,target,umis,genes" > $indir/nextgem_3p_satcurves.csv
for directory in $indir/F*; do
    sample=$(basename $directory)
    for target in ${targets[@]}; do
        results=$directory/$target/$sample/outs/per_sample_outs/$sample/metrics_summary.csv;
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
        echo "NextGEM_3P,$sample,$target,$umi,$genes" >> $indir/nextgem_3p_satcurves.csv
    done
done
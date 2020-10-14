main_phenotype_path=.csv # table with phenotype data from UK Biobank

output_phenoFile=.tsv # output table with selected phenotype data

ukb_cols_path=.txt # names of columns with phenotypes of interest

n="$(grep -c . $ukb_cols_path)"
readarray cols < $ukb_cols_path

sed -e 's/","/\t/g' -e 's/"//g' ${main_phenotype_path} | 
awk -F"\t"  -v n="$n" -v cols="${cols[*]}" -v OFS='\t' '{split(cols,list," "); for(i=1; i<=n; i++) {printf "%s%s", $list[i], OFS}; print "" }' > ${output_phenoFile}

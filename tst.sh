awk '{(if $1 == "rs10129287:106992524") print $0}' <(gzip -dc CHR_14.csv.gz)

awk -v FS=, '{if($1 == "rs10129287:106992524" && ) {print $0}}' <(gzip -dc CHR_14.csv.gz)
#!/usr/bin/env bash
set -euo pipefail
hm3=${1:?Usage: map.sh w_hm3.snplist output.tsv}
output=${2:?}
: "${PVAR_REF:?directory containing chr1.b38.pvar through chr22.b38.pvar}"
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
for chr in {1..22}; do
  awk -v chr="$chr" 'BEGIN{FS=OFS="\t"}
    FNR==NR{if(FNR>1)snp[$1]=1;next}
    $1 !~ /^#/ && ($3 in snp){print chr,$2,$3}' "$hm3" "$PVAR_REF/chr$chr.b38.pvar" >> "$tmp"
done
awk 'BEGIN{FS=OFS="\t";print "CHR","POS","SNP"}
  {key=$1":"$2;n[key]++;row[key]=$0}
  END{for(key in row)if(n[key]==1)print row[key]}' "$tmp" > "$output"

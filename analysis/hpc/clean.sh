#!/usr/bin/env bash
set -euo pipefail
format=${1:?Usage: clean.sh neale37|ea37|move37|move38|gcst38 input output}
input=${2:?}; output=${3:?}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/move.clean.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
lift_dir=${LIFT_DIR:-$PWD/liftover}
lift=$lift_dir/liftOver; chain=$lift_dir/hg19ToHg38.over.chain.gz
imp_dir=${PVAR_REF:-}; maf=.01
mkdir -p "$(dirname "$output")" "$lift_dir"
name=$(basename "$output" .gz)

download() {
  local url=$1 out=$2
  [[ -s $out ]] || curl -L --retry 3 -C - -o "$out" "$url"
}

need_liftover() {
  [[ -x $lift ]] || {
    download http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver "$lift"
    chmod +x "$lift"
  }
  [[ -s $chain ]] || download https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz "$chain"
}

lift37to38() {
  local in=$1 out=$2 name=$3 p
  p=$tmp/$name
  need_liftover

  gzip -cd "$in" | awk -v p="$p" '
    BEGIN { FS=OFS="\t" }
    NR==1 { print > p".header"; next }
    $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ {
      id=NR-1
      print "chr"$1,$2-1,$2,id > p".bed"
      print id,$0 > p".dat"
    }'

  "$lift" "$p.bed" "$chain" "$p.map" "$p.unmap" >/dev/null

  {
    cat "$p.header"
    awk '
      BEGIN { FS=OFS="\t" }
      FNR==NR {
        chr[$4]=$1; sub(/^chr/,"",chr[$4]); pos[$4]=$3
        next
      }
      ($1 in pos) {
        $2=chr[$1]; $3=pos[$1]
        if ($4 !~ /^rs/) $4=$2":"$3":"$6":"$5
        print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11
      }' "$p.map" "$p.dat"
  } | gzip -f > "$out"

  rm -f "$p".{header,bed,dat,map,unmap}
}

clean_neale37() {
  local in=$1 out=$2 name=$3 p
  p=$tmp/$name.hg19.gz
  gzip -cd "$in" | awk '
    BEGIN { FS=OFS="\t"; print "CHR","POS","SNP","EA","NEA","EAF","BETA","SE","P","N" }
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    $(h["low_confidence_variant"])!="true" &&
    $(h["n_complete_samples"])+0>0 &&
    $(h["beta"])!="NaN" && $(h["se"])!="NaN" && $(h["pval"])!="NaN" {
      split($1,a,":")
      chr=a[1]; pos=a[2]; ref=a[3]; alt=a[4]; n=$(h["n_complete_samples"])
      eaf=$(h["AC"])/(2*n)
      if (chr ~ /^[0-9]+$/ && chr>=1 && chr<=22)
        print chr,pos,$1,alt,ref,eaf,$(h["beta"]),$(h["se"]),$(h["pval"]),n
    }' | gzip -f > "$p"
  lift37to38 "$p" "$out" "$name"
  rm -f "$p"
}

clean_ea37() {
  local in=$1 out=$2 name=$3 p
  p=$tmp/$name.hg19.gz
  gzip -cd "$in" | awk '
    BEGIN { FS=OFS="\t"; print "CHR","POS","SNP","EA","NEA","EAF","BETA","SE","P","N" }
    NR==1 { next }
    $7!="NA" && $8!="NA" && $10!="NA" && $2 ~ /^[0-9]+$/ && $2>=1 && $2<=22 {
      print $2,$3,$1,$4,$5,$6,$7,$8,$10,765283
    }' | gzip -f > "$p"
  lift37to38 "$p" "$out" "$name"
  rm -f "$p"
}

clean_move37() {
  local in=$1 out=$2 name=$3 p
  p=$tmp/$name.hg19.gz
  gzip -cd "$in" | awk '
    BEGIN { FS=OFS="\t"; print "CHR","POS","SNP","EA","NEA","EAF","BETA","SE","P","N" }
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
      chr=$(h["CHR"]); eaf=$(h["EAF"])+0; maf=eaf<0.5 ? eaf : 1-eaf
      if (chr ~ /^[0-9]+$/ && chr>=1 && chr<=22 &&
          maf>0.01 && $(h["BETA"])!="NA" && $(h["SE"])!="NA" && $(h["P"])!="NA")
        print chr,$(h["POS"]),$(h["SNP"]),$(h["EA"]),$(h["NEA"]),$(h["EAF"]),$(h["BETA"]),$(h["SE"]),$(h["P"]),$(h["N"])
    }' | gzip -f > "$p"
  lift37to38 "$p" "$out" "$name"
  rm -f "$p"
}

clean_move38_pvar() {
  local in=$1 out=$2 name=$3 p chr dat
  p=$tmp/$name
  rm -f "$p".chr*.dat

  gzip -cd "$in" | awk -v p="$p" -v maf="$maf" '
    BEGIN { FS=OFS="\t" }
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
      chr=$(h["CHR"]); eaf=$(h["EAF"])+0; m=eaf<0.5 ? eaf : 1-eaf
      if (chr ~ /^[0-9]+$/ && chr>=1 && chr<=22 &&
          m>maf && $(h["BETA"])!="NA" && $(h["SE"])!="NA" && $(h["P"])!="NA")
        print $(h["SNP"]),$(h["EA"]),$(h["NEA"]),$(h["EAF"]),$(h["BETA"]),$(h["SE"]),$(h["P"]),$(h["N"]) >> p".chr"chr".dat"
    }'

  {
    echo -e "CHR\tPOS\tSNP\tEA\tNEA\tEAF\tBETA\tSE\tP\tN"
    for chr in {1..22}; do
      dat=$p.chr$chr.dat
      [[ -s $dat ]] || continue
      awk '
        BEGIN { FS=OFS="\t" }
        FNR==NR { rec[$1]=$0; next }
        $1 !~ /^#/ && ($3 in rec) { print $1,$2,rec[$3] }
      ' "$dat" "$imp_dir/chr$chr.b38.pvar"
    done
  } | gzip -f > "$out"

  rm -f "$p".chr*.dat
}

clean_gcst38() {
  local in=$1 out=$2
  gzip -cd "$in" | awk '
    BEGIN { FS=OFS="\t"; print "CHR","POS","SNP","EA","NEA","EAF","BETA","SE","P","N" }
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
      chr=$(h["chromosome"]); pos=$(h["base_pair_location"])
      snp=$(h["variant_id"]); if (snp=="" || snp==".") snp=$(h["rsid"])
      beta=$(h["beta"]); se=$(h["standard_error"]); p=$(h["p_value"]); n=$(h["n"])
      if (chr ~ /^[0-9]+$/ && chr>=1 && chr<=22 && beta!="." && se!="." && p!="." && n!=".")
        print chr,pos,snp,$(h["effect_allele"]),$(h["other_allele"]),$(h["effect_allele_frequency"]),beta,se,p,n
    }' | gzip -f > "$out"
}

case "$format" in
  neale37) clean_neale37 "$input" "$output" "$name" ;;
  ea37) clean_ea37 "$input" "$output" "$name" ;;
  move37) clean_move37 "$input" "$output" "$name" ;;
  move38) test -n "$imp_dir"; clean_move38_pvar "$input" "$output" "$name" ;;
  gcst38) clean_gcst38 "$input" "$output" ;;
  *) exit 2 ;;
esac
gzip -t "$output"

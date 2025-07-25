#!/usr/bin/env bash
min_depth_of_coverage=$1 ### default 10
heterozygous_perc=$2 ### default 0.2
output_prefix=$3 ### name of output vcf file
snps_list_file=$4 ### spia vcf file (in AnnotationFiles)
input_file=$5 ### .snps file output from PaCBAM
analysis_id=$6

(
  echo "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${analysis_id}";
  awk -v htperc="${heterozygous_perc}" -v min_depth_of_coverage="${min_depth_of_coverage}" \
  -v snps_list_file="${snps_list_file}" '
BEGIN {
  OFS          = "\t";
  s            = ":";
  htperc_top   = 1 - htperc;
  print_all    = 1;
  fixed_fields = "." OFS "." OFS "." OFS "GT";
}
{
  if (FILENAME == snps_list_file && $0 !~ /^(#|CHR)/) {
    ids[$1 s $2 s $3 s $4 s $5]   = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS fixed_fields OFS "./.";
    print_all = 0;
  } else {
    key  = $fn["chr"] s $fn["pos"] s $fn["rsid"] s $fn["ref"] s $fn["alt"];
    if (FNR == 1) {
      for (i = 1; i <= NF; i++) {
        fn[$i] = i;
      }
    } else if (print_all || (key in ids && !(key in snp_counts))) {
      snp_counts[key]++;
      af = $fn["af"];
      gt = "./.";
      if ($fn["cov"] >= min_depth_of_coverage) {
        if (af < htperc) {
          gt = "0/0";
        } else if (af > htperc_top) {
          gt = "1/1";
        } else {
          gt = "0/1";
        }
      }
      print $fn["chr"], $fn["pos"], $fn["rsid"], $fn["ref"], $fn["alt"], fixed_fields, gt;
    }
  }
}
END {
  if (!print_all) {
    for (id in ids) {
      if (!(id in snp_counts)) {
        print ids[id];
      }
    }
  }
}' "${snps_list_file}" "$input_file" | sort -k1,1V -k2,2n
) > "${output_prefix}.vcf"

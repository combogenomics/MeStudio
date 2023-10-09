#!/bin/bash

## installation_dir holds the entire path to the scripts and binaries
## You either install automatically, or manually.
## When installing manually make sure you also add the installation directory
## to your PATH [https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path]
## Linux script provided by Univeristy of Florence; The Florence Computational Biology Group;
## https://www.bio.unifi.it/vp-175-our-research.html

INSTALLATION_DIR=
replacer="$INSTALLATION_DIR/ms_replacR.py"
mscheck="$INSTALLATION_DIR/mscheck"
msmine="$INSTALLATION_DIR/msmine"
msfasta="$INSTALLATION_DIR/msfasta"
msmatch="$INSTALLATION_DIR/msmatch"
msx="$INSTALLATION_DIR/msx"
analyzer="$INSTALLATION_DIR/ms_analyzR.py"
circ="$INSTALLATION_DIR/ms_circ.R"

POSITIONAL_ARGS=()
print_help () {
  echo "Usage: mestudio -f <str> -anno <str> -smart <str> -mo <str> -o <str>"
  echo
  echo "-f    <str>           genomic sequence file"
  echo "-anno    <str>           genomic annotation file"
  echo "-smart   <str>           methylated base calls file"
  echo "-mo   <str>           newline delimited motifs list"
  echo "-o  <str>           output directory"
  echo
}

argc=0
while [[ $# -gt 0 ]]; do
  case $1 in
    -f)
      GENOMIC_FASTA=$(greadlink -m "$2")
      ((argc++))
      shift
      shift
      ;;
    -anno)
      GENOMIC_GFF=$(greadlink -m "$2")
      shift
      shift
      ((argc++))
      ;;
    -smart)
      METHYLATION_GFF=$(greadlink -m "$2")
      shift
      shift
      ((argc++))
      ;;
    -mo)
      MOTIFS_FILE=$(greadlink -m "$2")
      shift
      shift
      ((argc++))
      ;;
    -o)
      OUTPUT_DIR=$(greadlink -m "$2")
      shift
      shift
      ((argc++))
      ;;
    -*|--*)
      echo "Unknown option $1"
      print_help
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift # past argument
      ;;
  esac
done

## print help if script is called with no args
if [ "$argc" -eq 0 ]; then 
  print_help
  exit 1
fi

set -- "${POSITIONAL_ARGS[@]}"

python3.9 "$replacer" -o "$OUTPUT_DIR" -anno "$GENOMIC_GFF" -f "$GENOMIC_FASTA" -smart "$METHYLATION_GFF"

base_genomic_gff=$(basename "$GENOMIC_GFF")
base_methylation_gff=$(basename "$METHYLATION_GFF")
base_genomic_fasta=$(basename "$GENOMIC_FASTA")

cd "$OUTPUT_DIR"

GENOMIC_GFF=$(greadlink -m "$base_genomic_gff")
METHYLATION_GFF=$(greadlink -m "$base_methylation_gff")
GENOMIC_FASTA=$(greadlink -m "$base_genomic_fasta")
# GENOMIC_GFF=$(*"_anno.gff")
# METHYLATION_GFF=$(*"_smart.gff")
# GENOMIC_FASTA=$(*"_genomic.gff")

"$mscheck" -g "$GENOMIC_GFF" -f "$GENOMIC_FASTA" -m "$METHYLATION_GFF" -o mscore --mo "$MOTIFS_FILE" --type CDS
"$msmine" mscore/params.ms
"$msfasta" mscore/params.ms
"$msmatch" mscore/params.ms
"$msx" mscore/params.ms

cd mscore
# rm *.ms
# # mapfile -t motifs< <(sed 's/\r//g' "$MOTIFS_FILE" | sort | uniq)
for motif in *; do
    cd $motif
    python3 "$analyzer" -o "msa" -CDS "$motif"_CDS.gff -nCDS "$motif"_nCDS.gff -tIG "$motif"_true_intergenic.gff -US "$motif"_upstream.gff -anno $GENOMIC_GFF -s
    cd msa
    Rscript "$circ"
    cd ../../
done

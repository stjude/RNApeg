#!/bin/bash
set -e -o pipefail

main() {

  #################
  # Download data #
  #################

  echo ""
  echo "=== Setup ==="
  echo "  [*] Downloading input files..." 
  dx-download-all-inputs --parallel > /dev/null

  ################
  # Housekeeping #
  ################

  echo "  [*] Performing some housekeeping..."

  echo "   [-] Classpath: $CLASSPATH"
  echo "   [-] Perl5lib:  $PERL5LIB"

  echo "   [*] Setting up reference files.."
  local_reference_dir=/stjude/reference
  if [ ! -d $local_reference_dir ]
  then 
    mkdir -p $local_reference_dir
  fi
  species="Homo_sapiens"
  #ref_name="GRCh37-lite"
  echo "      - Species: $species"
  echo "      - Genome build: $ref_name"
  ref_dir=$local_reference_dir/$species/$ref_name
  if [ ! -d $ref_dir ] 
  then 
    mkdir -p $ref_dir
  fi

  # Download global reference data
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/FASTA
  dx download -o $ref_dir -r project-F5444K89PZxXjBqVJ3Pp79B4:/global/reference/$species/$ref_name/mRNA

  # Moving files where cloud configs expects them
  mv $bam_path /home/dnanexus/${bam_name}
  mv $bam_index_path /home/dnanexus/${bam_index_name}

  ##########
  # RNApeg #
  ##########

  echo "=== RNApeg ==="

  echo " [*] Running junction_extraction_wrapper.pl"
  mkdir output
  if [ ${ref_name} == "GRCh37-lite" ]
  then
    docker run -v /home/dnanexus:/data -v /stjude/reference:/reference -v /home/dnanexus/output:/results ghcr.io/stjude/rnapeg:RNAPEG_VERSION -b /data/${bam_name} -O /results -f /reference/Homo_sapiens/${ref_name}/FASTA/${ref_name}.fa -r /reference/Homo_sapiens/${ref_name}/mRNA/Combined/all_refFlats.txt -rg /reference/Homo_sapiens/${ref_name}/mRNA/Combined/all_refFlats.txt
  else
    docker run -v /home/dnanexus:/data -v /stjude/reference:/reference -v /home/dnanexus/output:/results ghcr.io/stjude/rnapeg:RNAPEG_VERSION -b /data/${bam_name} -O /results -f /reference/Homo_sapiens/${ref_name}/FASTA/${ref_name}.fa -r /reference/Homo_sapiens/${ref_name}/mRNA/RefSeq/refFlat-sharp.txt -rg /reference/Homo_sapiens/${ref_name}/mRNA/RefSeq/refFlat-sharp.txt
  fi

  junctions=$(dx upload /home/dnanexus/output/${bam_name}.junctions.tab --brief)
  dx-jobutil-add-output junctions "$junctions" --class=file
  junctions_shifted=$(dx upload /home/dnanexus/output/${bam_name}.junctions.tab.shifted.tab --brief)
  dx-jobutil-add-output junctions_shifted "$junctions_shifted" --class=file
  junctions_shifted_bed=$(dx upload /home/dnanexus/output/${bam_name}.junctions.tab.shifted.bed --brief)
  dx-jobutil-add-output junctions_shifted_bed "$junctions_shifted_bed" --class=file
  junctions_annotated=$(dx upload /home/dnanexus/output/${bam_name}.junctions.tab.shifted.tab.annotated.tab --brief)
  dx-jobutil-add-output junctions_annotated "$junctions_annotated" --class=file

}

#!/bin/bash
# Performs setup tasks for running RNApeg. 
#
# Accepts the following parameters: 
# $1 = target
# $2 = genome
# $3 = analysis configuration file
# $4 = run directory
# $5 = data directory
# $6 = output directory

# Show usage information if no parameters were sent
if [ "$#" == 0 ]; then about.sh $0; exit 1; fi

# Get parameters
ARGS=
while [ "$#" -gt 0 ]
do
  case "$1" in
  -*)             echo "Unrecognized option: $1" >&2; exit 1 ;;
  *)              ARGS="$ARGS $1"
  esac
  shift
done
set -- $ARGS
TARGET=$1
GENOME=$2
ANLS_CONFIG=$3
RUN_DIR=$4
DATA_DIR=$5
OUT_DIR=$6 

cp $ANLS_CONFIG $RUN_DIR/config.txt

# Source the 3 relevant config files in order.
echo
echo "Reading config files (you may see some error messages--these are OK):"
echo "* Application level"
echo "* Sequencing type level"
echo "..."
if [ `which_config.sh app rnapeg` ]
then
  . import_config.sh app rnapeg
fi 
if [ `which_config.sh target $TARGET` ]
then 
  . import_config.sh target $TARGET
fi
. import_config.sh genome $GENOME

. steplib.sh

set_step_script_dir $RUN_DIR 

init_step extract_junctions
cat > `get_step_make_cmds_script` <<EOF
#!/bin/bash
touch `get_step_cmds_file`
cat /dev/null > `get_step_cmds_file`

while read case_bam
 do
   bam="$DATA_DIR/\$case_bam/\$case_bam.bam"
   mkdir -p $OUT_DIR/\$case_bam
   echo "junction_extraction_wrapper.pl -bam \$bam -o $OUT_DIR/\$case_bam -genome $GENOME -now" >> `get_step_cmds_file`
 done < $RUN_DIR/config.txt
EOF
write_step_submit_script

init_step final_qa
cat > `get_step_qc_script` <<EOF
#!/bin/bash
# QC:
echo -n "" > final_qa.txt
anyfail=no
while read case_bam
do
  if ! qcquiet.sh `get_step_failed_qc_dir`/\$case_bam qc_rnapeg.sh $OUT_DIR/ \$case_bam
  then 
    anyfail=yes
    echo "FAIL \$case_bam" >> final_qa.txt
  else
    echo "PASS \$case_bam" >> final_qa.txt 
  fi
done < $RUN_DIR/config.txt
if [ "\$anyfail" == "yes" ]
then 
  echo There were QA failures
  echo "Exiting..."
  exit 1
fi
EOF



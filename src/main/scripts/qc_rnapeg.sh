#!/bin/bash
# QCs the RNApeg results for a run
#
# $1 = bam 
# $2 = data directory
source qclib.sh

DIR=$1 
BAM=$2

starttestcase RNApeg Using dir $DIR/$BAM

starttest DirExists
if [ -d "$DIR/$BAM" ]; then passtest; else aborttestcase; fi
echo $DIR/$BAM
cd $DIR/$BAM

starttest AnyFilesExist
if ls $DIR/$BAM/*.tab $DIR/$BAM/*.bed >&2; then passtest; else aborttestcase; fi

starttest AllFilesExist
if [ $(ls $DIR/$BAM/*.junctions.tab $DIR/$BAM/*.shifted.tab $DIR/$BAM/*.shifted.bed $DIR/$BAM/*.annotated.tab | wc -l) -ne '4' ]
then
	aborttestcase
else
	passtest
fi
 
starttest RNApegOutputSize 
cutoff=0
while read file
do 
  if [ "`filesize $file`" -lt $cutoff ]
  then 
    ls -l $file >&2
    failtestifactive Found at least one output that was too small
  fi
done< <(ls $DIR/$BAM/*.tab $DIR/$BAM/*.bed | sed '$d' )
passtestbydefault

summarize

#!/bin/bash

set -euo pipefail

CURRENT_DIR=$PWD
WORKING_DIR=$CURRENT_DIR/working_dir
OUTPUT_DIR=$CURRENT_DIR/output_dir
DATA_DIR=$CURRENT_DIR/data

APYTRAM=$PWD/../apytram.py

echo "START:"
date
START=$(date +%s)

mkdir -p $WORKING_DIR
cd $WORKING_DIR

$APYTRAM -d $WORKING_DIR/db/souris:MM,$WORKING_DIR/db/hamster:HM -dt paired:MM,paired:HM -out $OUTPUT_DIR/apytram -fq $DATA_DIR/rna_seq/fastq/Mesocricetus_auratus.datatest.fq:HM,$DATA_DIR/rna_seq/fastq/Mus_musculus.datatest.fq:MM -q $DATA_DIR/query/reference.nogap.fa:HH_Mus --out_by_species  -log apytram.log


files="$OUTPUT_DIR/apytram.HH_Mus.best.fasta
$OUTPUT_DIR/apytram.HH_Mus.HM.best.fasta
$OUTPUT_DIR/apytram.HH_Mus.MM.best.fasta
$OUTPUT_DIR/apytram.HH_Mus.fasta
$OUTPUT_DIR/apytram.HH_Mus.HM.fasta
$OUTPUT_DIR/apytram.HH_Mus.MM.fasta
"
for file in $files
do
if [ -s "$file" ]
then
    l="`wc -l < $file`"
    if [ "$l" -gt "1" ]
    then
        echo "$file has some data. -> OK"
    else
        echo "$file has no data."
        exit 1
    fi
else
    echo "$file is empty."
    exit 1
fi
done

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
echo "END:"
date


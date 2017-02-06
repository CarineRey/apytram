#!/bin/bash

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

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
echo "END:"
date


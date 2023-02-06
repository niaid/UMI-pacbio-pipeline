#!/bin/bash

forward=GAGCAGAAGACAGTGGCAATGA
forwardPrimer_RC=TCATTGCCACTGTCTTCTGCTC

reverse=CCCGCGTGGCCTCCTGAATTAT
reversePrimer_RC=ATAATTCAGGAGGCCACGCGGG

RTprimer=TATAATAAATCCCTTCCAGTCCCCCC
RTprimer_RC=GGGGGGACTGGAAGGGATTTATTATA

error=0.1
minLength=2800
maxLength=4000

referenceFile="reference-files/SHIV-AD8-EO_env.fasta"
dbFile="shiv_db/SHIV-AD8-EO_full-length.fasta"
term="SHIV"
umiLength=8

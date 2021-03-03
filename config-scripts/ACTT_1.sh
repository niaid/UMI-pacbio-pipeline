#!/bin/bash

## below is the information for ACTT human cells covid-19 for SEM_protein region

forward=TTCGCATGGTGGACAGCCTTTGTT
forwardPrimer_RC=AACAAAGGCTGTCCACCATGCGAA

reverse=CCGCTCCGTCCGACGACTCACTATACCCGCGTGGCCTCCTGAATTAT
reversePrimer_RC=ATAATTCAGGAGGCCACGCGGGTATAGTGAGTCGTCGGACGGAGCGG

RTprimer=TCTCCATTGGTTGCTCTTCATCT
RTprimer_RC=AGATGAAGAGCAACCAATGGAGA


error=0.30
minLength=4500
maxLength=6500

referenceFile="/hpcdata/vrc_vpds/data/pacbio/pacbio-pipeline/reference-files/CoV_Reference_SEM_protein.fa"

#!/bin/bash

## below is the information for invitro monkey cells covid-19 fo SEM_protein region

forward=TTCGCATGGTGGACAGCCTTTGTT        
forwardPrimer_RC=AACAAAGGCTGTCCACCATGCGAA

reverse=CCGCTCCGTCCGACGACTCACTATA        
reversePrimer_RC=TATAGTGAGTCGTCGGACGGAGCGG

RTprimer=TCTCCATTGGTTGCTCTTCATCT
RTprimer_RC=AGATGAAGAGCAACCAATGGAGA

error=0.30
minLength=4500
maxLength=6500

referenceFile="/hpcdata/vrc_vpds/data/pacbio/pacbio-pipeline/reference-files/CoV_Reference_SEM_protein.fa"

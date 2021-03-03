#!/bin/bash

## below is the information for invitro monkey cells covid-19 fo S_protein region

forward=TTCGCATGGTGGACAGCCTTTGTT        
forwardPrimer_RC=AACAAAGGCTGTCCACCATGCGAA

reverse=CCGCTCCGTCCGACGACTCACTATA        
reversePrimer_RC=TATAGTGAGTCGTCGGACGGAGCGG

RTprimer=CGTTGCAGTAGCGCGAACAA
RTprimer_RC=TTGTTCGCGCTACTGCAACG

error=0.30
minLength=2800
maxLength=5000

referenceFile="reference-files/CoV_Reference_S_protein.fa"

#!/bin/bash

forward=TTCGCATGGTGGACAGCCTTTGTT        
forwardPrimer_RC=AACAAAGGCTGTCCACCATGCGAA

reverse=CCGCTCCGTCCGACGACTCACTATA        
reversePrimer_RC=TATAGTGAGTCGTCGGACGGAGCGG

RTprimer=CGTTGCAGTAGCGCGAACAA
RTprimer_RC=TTGTTCGCGCTACTGCAACG

error=0.10
minLength=2800
maxLength=5000

referenceFile="reference-files/CoV_Reference_S_protein.fa"
dbFile="cov_db/COVdb.fa"
term="coronavirus"
umiLength=8

# uodeer
Draft transcriptome and analysis of deer for protein identification and quantification

**NOTE: It is assumed that the reader of this has a basic understanding of command line tools, can navigate around and install packages.**

>The ultimate goal of this is to generate a draft transcriptome. This can then have putative ORFs translated and give rise to a peptide identification database for deer. This will then be compared to previously obtained mass spec data to give us hopefully an increased number of proteins identified.

>**Firstly evaluating the raw reads**

```bash
java -jar fastqc.jar ~/path/to/raw/seqs/*.fq.gz 

```

>This will generate fastQC checks for all sequences that you can go through to check and make sure the original data is actually good quality.


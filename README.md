# uodeer
Draft transcriptome and analysis of deer for protein identification and quantification

**NOTE: It is assumed that the reader of this has a basic understanding of command line tools, can navigate around and install packages.**

>The ultimate goal of this is to generate a draft transcriptome. This can then have putative ORFs translated and give rise to a peptide identification database for deer. This will then be compared to previously obtained mass spec data to give us hopefully an increased number of proteins identified.

>**Firstly evaluating the raw reads**

```bash
java -jar fastqc.jar ~/path/to/raw/seqs/*.fq.gz 

```

>This will generate fastQC checks for all sequences that you can go through to check and make sure the original data is actually good quality.

>**Trimming the raw reads using Trimmomatic**

```bash
java -jar trimmomatic-0.35.jar PE -threads 20 -phred33 \
 sample_read1.fastq.gz sample_read2.fastq.gz \
 output1_forward_paired.fq.gz output1_forward_unpaired.fq.gz \
 output1_reverse_paired.fq.gz output1_reverse_unpaired.fq.gz \
 ILLUMINACLIP:~/path/to/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

```

>This is quite lenient trimming, but i have done so to ensure we keep as much RNA as possible. * Note copying and pasting this code wont make it work, you will need to link directly to the correct place for the trimmomatic.jar and also link to the adapters. For bioc this is /usr/local/bin/Trimmomatic-0.35/trimmomatic.jar and /usr/local/bin/Trimmomatic-0/35/adapterss/TruSeq3-PE.fa respectivley.

>**Post trim QC**

>Just run fastQC on the trimmed paired forward and reverse files, to confirm the removal of adapter sequence and low qual reads.

Once sample pre-processing is complete, trinity assemblies can be carried out.

>**Trinity assemblies**

```bash
trinity --CPU 20 --seqType fq --max_memory 50G --left output1_forward_paired.fq.gz \
--right output1_reverse_paired.fq.gz --output trinity_output1

```

>Trinity was run on all samples individually and then on all the samples at once. To add more than one sample to the assembly just put the reads seperated by a comma (no space) eg; --left output1_forward_paried.fq.gz,output2_forward_paried.fq.gz . . .






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

> **Post trim QC**

>Just run fastQC on the trimmed paired forward and reverse files, to confirm the removal of adapter sequence and low qual reads.

Once sample pre-processing is complete, trinity assemblies can be carried out.

> **Trinity assemblies**

```bash
trinity --CPU 20 --seqType fq --max_memory 50G --left output1_forward_paired.fq.gz \
--right output1_reverse_paired.fq.gz --output trinity_output1

```

>Trinity was run on all samples individually and then on all the samples at once. To add more than one sample to the assembly just put the reads seperated by a comma (no space) eg; --left output1_forward_paried.fq.gz,output2_forward_paried.fq.gz . . .

From the Trinity assemblies you get a .fasta file that contains A LOT of different predicted transcripts. It is not fesible to process all of these. The easiest way to deal with this is to re-map the trimmed files back onto the .fasta file and do raw read counting. In theory the most expressed transcripts (those with the highest counts) are also going to be some of the proteins we care about (although this in not always the case, for the first round of analysis its much easier to simplify it this way).

In my experiences the Trinity quant pipeline is not the greatest. Essentially the way it parses your mapped file is the wrong way for most programs that generate counts (they sort the BAM file relative to the reference.fasta, whereas these programs need a name sorted bam file). So using Trinity for the first step of the mapping is alright to do, but the actual read counting will probably fall over.

> **Trinity mapping/counting**

>Essentially this will map the reads to the .fasta file generated from Trinity. More importantly it will also map the reads to everywhere they match (allows unlimited non-unique mapping). This is important as it allows the count software to pick which transcripts are most likely, without being biased by half your reads mapping to a truncated transcript. This is what a lot of people used to qunatify their transcripts.

```bash
/Trinity/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq \
--left output1_forward_paired.fq.gz --right output1_reverse_paired.fq.gz --est_method salmon \
--aln_method bowtie2 --trinity_mode --prep_reference --output_dir eXpress_output

```

>This might error out during the read counting of salmon. My solution to this - Just use salmon by itself

```bash
salmon -quant -i indexed_Trinity.fasta -l IU -1 output1_forward_paired.fq.gz \
-2 output1_reverse_paired.fq.gz -o salmon_quant_output1 -p 20

```

>This SHOULD work, if trinity fell over.

>From here the most abundantly expressed transcripts (in prelim analysis it was top 1000) were listed from the "quant.sf" file from salmon. However, the sequences need to be extracted from the original Trinity.fasta assembly so they can be translated. This is a really handy way of doing it - Make all of the transcripts you want into a list in a text file (make sure it has unix coding not microsoft, otherwise this wont work, less say this is call top1000list.txt).

```bash
cat top1000list.txt | xargs samtools faidx Trinity.fasta | cat > top1000transcripts.fasta

```

>This will use the fasta header identifier to rapidly pull them out and put them all into a file of their own (very handy)

>The next step is to take these and turn them into logical putative proteins. I did this initially in 2 different ways. 1) simply translating anything that had an AUG and an inframe stop codon that was over 50 residues away (brute force method). and 2) A measured approach using hmm and blast matches.

> **Identifying putative transcripts using TRANSdecoder

```bash
transDecoder.LongOrfs -t top1000transcripts.fasta -m 50

```

>This will output quite a lot of transcripts, although we want some that have certainty, to do this i incorperated BLASTp searches and hmm comparisons to determine proteins that are more likely. In this case all BLASTp was done in comparison to Bos_taurus only and hmm was compared to all of Pfam.

```bash
blastp -query top1000transcripts/longest_orfs.pep  -db bos_taurus  -max_target_seqs 1 \
-outfmt 6 -evalue 1e-5 -num_threads 20 > blastp.outfmt6

hmmscan --cpu 20 --domtblout pfam.domtblout Pfam-A.hmm top1000transcripts/longest_orfs.pep

TransDecoder.Predict -t Trinity.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

```

>The output of longest_orfs.pep can just be changed into fasta if you remove spaces from the fasta > headers (some programs really dont like spaces).













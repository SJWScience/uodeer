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

From the Trinity assemblies you get a .fasta file that contains A LOT of different predicted transcripts. It is not fesible to process all of these. The easiest way to deal with this is to re-map the trimmed files back onto the .fasta file and do raw read counting. In theory the most expressed transcripts (those with the highest counts) are also going to be some of the proteins we care about (although this in not always the case, for the first round of analysis its much easier to simplify it this way).

In my experiences the Trinity quant pipeline is not the greatest. Essentially the way it parses your mapped file is the wrong way for most programs that generate counts (they sort the BAM file relative to the reference.fasta, whereas these programs need a name sorted bam file). So using Trinity for the first step of the mapping is alright to do, but the actual read counting will probably fall over.

>** Trinity mapping/counting**

>Essentially this will map the reads to the .fasta file generated from Trinity. More importantly it will also map the reads to everywhere they match (allows unlimited non-unique mapping). This is important as it allows the count software to pick which transcripts are most likely, without being biased by half your reads mapping to a truncated transcript.

```bash
/Trinity/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq \
--left output1_forward_paired.fq.gz --right output1_reverse_paired.fq.gz --est_method eXpress \
--aln_method bowtie2 --trinity_mode --prep_reference --output_dir eXpress_output

```

>This will probably error out during the read counting of eXpress with some error about alignments not being properly sorted. My solution to this - Run eXpress by itself, after sorting the bowtie2.bam file that gets output from the above code.

```bash

samtools sort -n reads.bam > reads.sorted.bam

express Trinity.fasta reads.sorted.bam -o eXpress_output1

```

>This SHOULD work, if trinity fell over. Like i mentioned above, i think that trinity script either doesnt sort the .bam file or sorts it via position not by name (samtools sort -n) like eXpress required.









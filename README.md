# JaRMS

A java port of the [CNVnator](http://www.ncbi.nlm.nih.gov/pubmed/21324876) Mean-shift algorithm for predicting CNVs. JaRMS stands for **Ja**va **R**ead depth **M**ean **S**hift.

## Installation

JaRMS requires Java version 1.8 to be installed on your system. 

In order to install, just download the git repository:

```bash
git clone https://github.com/njdbickhart/JaRMS
```

You can currently invoke JaRMS as a "jar" file within the repository as we work to implement all features. In order to run JaRMS this way, you simply have to type the following java command:

```bash
java -jar ./JaRMS/store/JaRMS.jar
```

You should see a "help" message if your java executable is installed correctly. 

## Running the program

As we develop JaRMS, there are some developer options currently enabled to assist with the interpretation of data. For the purposes of using the software for research purposes, we only support the use of the "call" mode at this time. In order to retrieve the help menu for "call" mode, just invoke the program on the command line:

```bash
java -jar ./JaRMS/store/JaRMS.jar call
JaRMS call mode
Usage: java -jar JaRMS.jar cluster [-i bamfile -f fasta file -o output prefix] (option: -t number of threads)
        -i      One or more BWA-processed bam files for processing
        -f      The reference genome fasta file that was used during alignment of the bam file
        -o      Output file prefix and directory
        -t      Number of threads to use [optional: use one thread]
        -w      Use this window size [optional: determine from BAM read depth]
```

Here is an explanation of the inputs to the program:

| Argument | Required? | Description |
| :--- | :--- | :--- |
| -i | yes | A BWA-aligned BAM file containing the sequence reads used for calling CNVs |
| -f | yes | The same reference fasta file used to align reads in the BAM file in "-i" |
| -o | yes | The output directory and file name prefix for all subsequent output files |
| -t | no | The number of threads to use (currently disabled! Only uses one thread) |
| -w | no | The number of base pairs used to generate non-overlapping windows for CNV calling |

Most arguments are self-explanatory; however, the window size is a crucial determinant for the success of your CNV calling. JaRMS will attempt to automatically estimate whether or not a window size of 500 bp or 100 bp is appropriate for your data. To override this estimate, specify a different value using "-w." 

## Interpreting the output

If JaRMS completes successfully, you should have two files:

* A BED file containing predicted CNV intervals and statistics to give you the likelihood of the CNV call
* A "levels" file containing condensed copy number estimates for segments of the genome for that individual

#### The BED file with CNV calls

Here is a description of the information contained in the CNV call file. The columns are listed in order of their appearance from left to right:

1. Chromosome/scaffold/contig
2. Start coordinate (bp)
3. End coordinate (bp)
4. Normalized read depth (with "1" as the base line read depth for the entire sample)
5. The p value from a T test of the average read depth in the region compared to genome variability estimates (assumed that only 10% of the genome is copy number variable)
6. The probability that the average read depth values for the region are within the tails of a normal distribution of read depth values across the entire genome


#### The LEVELS file

Here is a description of the information contained in the Levels file. This file contains normalized read depth information intervals that were condensed by the Mean Shift algorithm. The columns are listed in order of their appearance from left to right:

1. Chromosome/scaffold/contig
2. Start coordinate (bp)
3. End coordinate (bp)
4. Normalized read depth for the region (with "1" as the base line read depth for the entire sample)
# sci-L3-target-seq primer designer

## Introduction


## Prerequisites
python=3.6  
pandas=1.0.1  
biopython  
bedtools=2.29.1  
primer3 (libprimer3 release 2.5.0)  
ispcr  

or create and activate environment in scripts directory  
>> $ conda env update -f ./scripts/environment.yaml  
>> $ conda activate sci-L3-target-seq-primers

## Usage
### Prepare input files
#### 1. Create a BED6 file containing the regions to be targetted.
        Tab-delimeted 6 columns with .bed extension:
        [chromosome, region start, region end, unique gene/exon name (without spaces), score(default 1 for all), strand(+ or -)]
        Example:
            chr9	133589332	133589842	ABL1_Ex1	1	+
            chr9	133729450	133729624	ABL1_Ex2	1	+
#### 2. Create a CSV file containing a list of barcodes for the universal primer
#### 3. Create a CSV file containing a list of sample specific indices


### Run
    $ python ./scripts/run.py -i <path to bed file> -g <path to fasta file of reference genome> -b <path to .CSV file of barcodes> -indices <path to .CSV file with sample specific indices> [options]
    
### Optional Arguuments

    -h, --help                  show this help message and exit
    -i, --input                 Enter the path to the BED or FASTA file
    -g, --genome                Enter the path to the reference genome
    --P5, -P5                   P5 Adapter sequence(Default=AATGATACGGCGACCACCGAGA)
    --P7, -P7                   P7 Adapter sequence(Default=CAAGCAGAAGACGGCATACGAGAT)
    --read1, -r1                Read1 sequence(Default=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)
    --read2, -r2                Read2 sequence(Default=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)
    --SSSPR, -ssspr             SSS Primer Region(Default=GGGATGCAGCTCGCTCCTG)
    --indices, -indices         Enter the path to the CSV file containing sample specific indices
    --barcodes, -b              Enter the path to the CSV file containing Barcodes

# STAMPS_Tutorial2017

## Getting started

This tutorial describes a complete walkthrough for CONCOCT2 and DESMAN on the STAMPS servers using 4 threads. Some of of the pre-processing steps are time consuming and not that informative. These have been pre-run and you will just copy across the output.
We will begin by constructing a new working directory in our home directories:

```bash
cd ~
mkdir CDTutorial
cd CDTutorial
```

The tutorial data consists of 32 synthetic samples of 500,000 2X150 bp paired end reads taken from 8 species and 15 genomes. Copy across the simulation config file just to have a look at it with the less command:

```bash
cp /class/stamps-shared/CDTutorial/genome_coverage.tsv .
less genome_coverage.tsv
```

This has format "Sample\tSpecies\Strain\Coverage\Relative frequency". We will now use CONCOCT/DESMAN to resolve these *de novo*.

The simulated reads and genomes are in the tar archives /class/stamps-shared/CDTutorial/Reads.tar.gz and /class/stamps-shared/CDTutorial/Genomes.tar.gz respectively. Please **do not** run the commands below. They indicate what you would do if you wanted to run the preprocessing steps yourself. 

```
cp /class/stamps-shared/CDTutorial/Reads.tar.gz
tar -xvzf Reads.tar.gz
cp /class/stamps-shared/CDTutorial/Genomes.tar.gz
tar -xvzf Genomes.tar.gz
```
**do not run the above**


## Running CONCOCT2 and DESMAN on the 32 sample mock

We need to load up the concoct and desman modules:

```basg
module load concoct
source /class/stamps-software/concoct_speedup_mp/bin/activate
module load desman
source /class/stamps-software/desman/bin/activate
```

We now describe how to perform a complete analysis binning and resolving strains on 
this synthetic community. Some of these steps are time consuming so we recommend 
that you copy across the output instead. This assumes that DESMAN and CONCOCT are installed and their paths set to the variables DESMAN and CONCOCT respectively e.g. (changing paths to your system set-up):

```bash
export CONCOCT=/automounts/classfs/classfs/stamps-software/concoct_speedup_mp/
export DESMAN=/automounts/classfs/classfs/stamps-software/desman/
export CDSCRIPTS=/class/stamps-shared/CDTutorial/scripts
```

We will also create a new variable pointing to our current working dir for all this analysis:

```bash
export METASIMWD=~/CDTutorial
```

The first step in the analysis is to assemble the reads. 

### Assembly

I previously assembled the reads using MEGAHIT 1.1.1 and default parameters:

```
module load megahit/1.0.6
ls Reads/*r1*gz | tr "\n" "," | sed 's/,$//' > r1.csv
ls Reads/*r2*gz | tr "\n" "," | sed 's/,$//' > r2.csv
megahit -1 $(<r1.csv) -2 $(<r2.csv) -t 4 -o Assembly 
```

This would take approximately 20 minutes. Please __do not run__ the above now. Instead copy across results that I ran earlier:

```bash
cp /class/stamps-shared/CDTutorial/Assembly.tar.gz .
tar -xvzf Assembly.tar.gz
```  

Evaluate the assembly quality with the script provided:
```
$CDSCRIPTS/contig-stats.pl < Assembly/final.contigs.fa 
```

Results in output:
```
sequence #: 13863	total length: 31472365	max length: 1015941	N50: 9774	N90: 706
```

The next step in CONCOCT is to cut up contigs and map the reads back onto them. Again 
__do not run__ any of the following:


```
cd Assembly
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m final.contigs.fa > final_contigs_c10K.fa
bwa index final_contigs_c10K.fa
cd ..

mkdir Map

for file in Reads/*.r1.fq.gz
do

    stub=${file%.r1.fq.gz}

    base=${stub##*/}

    echo $base

    bwa mem -t 8 Assembly/final_contigs_c10K.fa $file ${stub}.r2.fq.gz > Map/$base.sam

done

python $DESMAN/scripts/Lengths.py -i Assembly/final_contigs_c10K.fa | tr " " "\t" > Assembly/Lengths.txt

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
    (samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam; samtools sort -m 1000000000 ${stub}.mapped.bam ${stub}.mapped.sorted; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt)
done
```

Collate coverages together (__do not run__):

```
for i in Map/*_cov.txt 
do 
   echo $i
   stub=${i%_cov.txt}
   stub=${stub#Map\/}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > Map/${stub}_cov.csv
done

```

Instead copy and extract out the prepared directory:
```
cp /class/stamps-shared/CDTutorial/Map.tar.gz .
tar -xvzf Map.tar.gz
```

### Contig binning

We are going to use a more efficient version of CONCOCT, CONCOCT2. Now we can run CONCOCT:

```
$DESMAN/scripts/Collate.pl Map > Coverage.csv

mkdir Concoct

mv Coverage.csv Concoct

cd Concoct

tr "," "\t" < Coverage.csv > Coverage.tsv

concoct --coverage_file Coverage.tsv --composition_file ../Assembly/final_contigs_c10K.fa -t 4

```

### Determine number of complete genomes

Find genes using prodigal:
```
    cd ..
    
    mkdir Annotate

    cd Annotate/

    python $DESMAN/scripts/LengthFilter.py ../Assembly/final_contigs_c10K.fa -m 1000 >     final_contigs_gt1000_c10K.fa

    prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d     final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff 
```

Assign COGs change the -c flag which sets number of parallel processes appropriately:
```
    export COGSDB_DIR=/class/stamps-shared/CDTutorial/rpsblast_cog_db/
    module load parallel/201702222
    $CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 4 -r 1
```

We will also write out lists of cogs and genes in each contig. These will be useful later:
```
python $DESMAN/scripts/ExtractCogs.py -b final_contigs_gt1000_c10K.out --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.cogs
python $DESMAN/scripts/ExtractGenes.py -g final_contigs_gt1000_c10K.gff > final_contigs_gt1000_c10K.genes
```

These are just simple comma separated lists, have a look at them:

```
head final_contigs_gt1000_c10K.cogs 
head final_contigs_gt1000_c10K.genes
```

Format is:

* gene_id,cog_id,startpos,endpos,cogstart,cogend,strand
* gene_id,contig_id,startpos,endpos,strand

We will also need contig lengths later on:
```
python $DESMAN/scripts/Lengths.py -i final_contigs_gt1000_c10K.fa > final_contigs_gt1000_c10K.len
```

We are going to determine the number of complete and pure genomes using single-core gene frequencies. First we 
calculate scg frequencies on the CONCOCT clusters:

```
cd ../Concoct
python $CONCOCT/scripts/COG_table.py -b ../Annotate/final_contigs_gt1000_c10K.out  -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c clustering_gt1000.csv  --cdd_cog_file $CONCOCT/scgs/cdd_to_cog.tsv > clustering_gt1000_scg.tsv
```

Then we visualise in R:
```
$CONCOCT/scripts/COGPlot.R -s clustering_gt1000_scg.tsv -o clustering_gt1000_scg.pdf
```

We should see 8 nearly complete bins and some fragmentary ones:

![SCG Results](Figures/clustering_gt1000_scg.pdf)

We want to store a list of 75% complete and pure for DESMAN analysis:
```
python $SCRIPTS/CompleteClusters.py clustering_gt1000_scg.tsv > Cluster75.txt
```

### How to calculate cluster coverages

```
python $EXAMPLE/scripts/ClusterMeanCov.py Coverage.csv clustering_gt1000.csv ../Assembly/final_contigs_c10K.fa > Cluster_Cov.csv
```

## Run DESMAN pipeline to resolve strains in each high quality bin


First we separate out the contigs, cogs, and genes into the separate bins:
```
cd ..
mkdir Split
cd Split
$DESMAN/scripts/SplitClusters.pl ../Annotate/final_contigs_gt1000_c10K.fa ../Concoct/clustering_gt1000.csv
$METASIMPATH/scripts/SplitCOGs.pl ../Annotate/final_contigs_gt1000_c10K.cogs ../Concoct/clustering_gt1000.csv
$METASIMPATH/scripts/SplitGenes.pl ../Annotate/final_contigs_gt1000_c10K.genes ../Concoct/clustering_gt1000.csv
cd ..
```

Now we can split up the bam files by each cluster in turn:
```
mkdir SplitBam

while read -r cluster 
do
    grep ">" Split/${cluster}/${cluster}.fa | sed 's/>//g' > Split/${cluster}/${cluster}_contigs.txt
    ./STAMPS_Tutorial2017/scripts/AddLengths.pl Annotate/final_contigs_gt1000_c10K.len < Split/${cluster}/${cluster}_contigs.txt > Split/${cluster}/${cluster}_contigs.tsv
    mkdir SplitBam/${cluster}

    for bamfile in Map/*.mapped.sorted.bam
    do
        stub=${bamfile#Map\/}
        stub=${stub%.mapped.sorted.bam}
        
        samtools view -bhL Split/${cluster}/${cluster}_contigs.tsv $bamfile > SplitBam/${cluster}/${stub}_Filter.bam&

    done 
    wait    
done < Concoct/Cluster75.txt 
```

and use a third party program bam-readcount to get base frequencies at each position on each contig:

```
while read line
do
    mag=$line

    echo $mag

    cd SplitBam
    cd ${mag}

    cp ../../Split/${mag}/${mag}_contigs.tsv ${mag}_contigs.tsv
    samtools faidx ../../Split/${mag}/${mag}.fa

    echo "${mag}_contigs.tsv"
    mkdir ReadcountFilter
    for bamfile in *_Filter.bam
    do
        stub=${bamfile%_Filter.bam}
            
        echo $stub

        (samtools index $bamfile; bam-readcount -w 1 -q 20 -l "${mag}_contigs.tsv" -f ../../Split/${mag}/${mag}.fa $bamfile 2> ReadcountFilter/${stub}.err > ReadcountFilter/${stub}.cnt)&
    
     done
     cd ..
     cd ..
     wait
    
done < Concoct/Cluster75.txt

``` 

Then we want to select core cogs from each cluster

```
while read -r cluster 
do
    echo $cluster
    $DESMAN/scripts/SelectContigsPos.pl STAMPS_Tutorial2017/data/cogs.txt < Split/${cluster}/${cluster}.cog > Split/${cluster}/${cluster}_core.cogs
done < Concoct/Cluster75.txt
```

Then we can get the base counts on these core cogs:

```
#!/bin/bash

mkdir Variants
while read -r cluster 
do
    echo $cluster
    cd ./SplitBam/${cluster}/ReadcountFilter
    gzip *cnt
    cd ../../..
    python $DESMAN/scripts/ExtractCountFreqGenes.py Split/${cluster}/${cluster}_core.cogs ./SplitBam/${cluster}/ReadcountFilter --output_file Variants/${cluster}_scg.freq > Variants/${cluster}log.txt&

done < Concoct/Cluster75.txt 
``` 



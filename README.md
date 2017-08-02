# STAMPS_Tutorial2017

This tutorial provides a complete walkthrough for CONCOCT2 and DESMAN on the STAMPS servers.
First download the reference genomes and simulation coverages used for validation of the results:

```
wget https://stamps2017.s3.climb.ac.uk/Genomes.tar.gz
tar -xvzf Genomes.tar.gz
wget https://stamps2017.s3.climb.ac.uk/coverage.tsv
wget https://stamps2017.s3.climb.ac.uk/Genomes.tar.gz
```

Now get the actual metagenome reads there are 32 samples:
```
mkdir Reads
cd Reads
```

Now enter this for loop as one command:
```
for n in `seq 1 32 `
do
    wget https://stamps2017.s3.climb.ac.uk/Reads.${n}.r1.fq.gz
    wget https://stamps2017.s3.climb.ac.uk/Reads.${n}.r2.fq.gz
done
cd..
```

## Running CONCOCT2 and DESMAN on the complex mock

We now describe how to perform a complete analysis binning and resolving strains on 
this synthetic community. Some of these steps are time consuming so we provide the option to download the output instead. This assumes that DESMAN and CONCOCT are installed and their paths 
set to the variables DESMAN and CONCOCT respectively e.g. (changing paths to your 
system set-up):

```
export CONCOCT=/mnt/gpfs/chris/repos/CONCOCT/
export DESMAN=/mnt/gpfs/chris/repos/DESMAN/
```

We will also create a new variable pointing to our current working dir for all this analysis:
```
export METASIMPATHWD=$METASIMPATH/ComplexStrainSim/Strains/Simulation
```

The first step in the analysis is to assemble the reads. 

### Assembly

We assembled the reads using MEGAHIT 1.1.1 and default parameters:
```
ls Reads/*r1*gz | tr "\n" "," | sed 's/,$//' > r1.csv
ls Reads/*r2*gz | tr "\n" "," | sed 's/,$//' > r2.csv
megahit -1 $(<r1.csv) -2 $(<r2.csv) -t 96 -o Assembly > megahit1.out
```

Then cut up contigs and index for BWA:

```
cd Assembly
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m final.contigs.fa > final_contigs_c10K.fa
bwa index final_contigs_c10K.fa
cd ..
```



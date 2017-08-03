# STAMPS_Tutorial2017

This tutorial provides a complete walkthrough for CONCOCT2 and DESMAN on the STAMPS servers using 
8 threads.
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
megahit -1 $(<r1.csv) -2 $(<r2.csv) -t 8 -o Assembly > megahit1.out
```

Evaluate the assembly quality with the script provided:
```
$SCRIPTPATH/contig-stats.pl < Assembly/final.contigs.fa 
```

Results similar to...
sequence #: 26972	total length: 70633244	max length: 1091259	N50: 13151	N90: 830
sequence #: 8678	total length: 31484306	max length: 1071063	N50: 15745	N90: 1406

Then cut up contigs and index for BWA:

```
cd Assembly
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m final.contigs.fa > final_contigs_c10K.fa
bwa index final_contigs_c10K.fa
cd ..
```

If you do not have time for the assembly simply download:
```
wget https://stamps2017.s3.climb.ac.uk/Assembly.tar.gz
tar -xvzf Assembly.tar.gz
```

```
mkdir Map

for file in Reads/*.r1.fq.gz
do

    stub=${file%.r1.fq.gz}

    base=${stub##*/}

    echo $base

    bwa mem -t 8 Assembly/final_contigs_c10K.fa $file ${stub}.r2.fq.gz > Map/$base.sam
done
```

### Generate contig coverages

And calculate coverages:
```
python $DESMAN/scripts/LengthsT.py -i Assembly/final_contigs_c10K.fa | tr " " "\t" > Assembly/Lengths.txt

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
    (samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam; samtools sort -m 1000000000 ${stub}.mapped.bam -o ${stub}.mapped.sorted.bam; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt)
done
```

for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub  
     bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g Assembly/Lengths.txt > ${stub}_cov.txt
done

### Contig binning

We are going to use a more efficient version of CONCOCT, CONCOCT2. Now we can run CONCOCT:
```

    mkdir Concoct

    mv Coverage.csv Concoct

    cd Concoct

    tr "," "\t" < Coverage.csv > Coverage.tsv

    concoct --coverage_file Coverage.tsv --composition_file ../Assembly/final_contigs_c10K.fa -t 8 > concoct.out

```

### Determine number of complete genomes

Find genes using prodigal:
```
    cd ..
    
    mkdir Annotate

    cd Annotate/

    python $DESMAN/scripts/LengthFilter.py ../Assembly/final_contigs_c10K.fa -m 1000 >     final_contigs_gt1000_c10K.fa

    prodigal -i final_contigs_gt1000_c10K.fa -a final_contigs_gt1000_c10K.faa -d     final_contigs_gt1000_c10K.fna  -f gff -p meta -o final_contigs_gt1000_c10K.gff > p.out
```

Assign COGs change the -c flag which sets number of parallel processes appropriately:
```
    export COGSDB_DIR=/mydatabase_path/rpsblast_db
    $CONCOCT/scripts/RPSBLAST.sh -f final_contigs_gt1000_c10K.faa -p -c 8 -r 1
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
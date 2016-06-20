#Canadian Bioinformatics Workshop 2016
#Marker Gene Integrated Assignment

#Make the script noisy
set -v

#Setup some parameters for the script
dataLocation="~/CourseData/metagenomics/integrated_assignment_day1/" # this is where your fastq files should be
workingDir="~/workspace/assignment1" # this is where all the files generated in this script will be saved
ncores=4 # how many cores do we have to work with

cd $workingDir
workingDir=$PWD #get absolute path

#First, we will link our sequence and reference data into our workspace
mkdir sequence_files
mkdir reference_data
mkdir scripts
ln -s $dataLocation/*.fastq.gz sequence_files/
ln -s $dataLocation/nwcs_arct_metadata.tsv .
ln -s $dataLocation/97* reference_data/
ln -s $dataLocation/*.py scripts

#Prep database for SortMeRNA using: cd reference_data; /usr/local/sortmerna-2.1/indexdb_rna --ref 97_otus.fasta
sortmernaDB=$workingDir"/reference_data/97_otus"

#Mark the scripts as executable
chmod u+x scripts/*

#Check the read files using FastQC - usually this is useful to check for length and low quality base calls
mkdir fastqc_out
mkdir fastqc_out/raw
mkdir fastqc_out/raw/individuals
fastqc -q -t $ncores sequence_files/*.fastq.gz -o fastqc_out/raw/individuals
mkdir fastqc_out/raw/combined
#May not need to look at all the individual files if the results from all the files combined looks *perfect*. Best to always look at your data, however.
gunzip --to-stdout sequence_files/*.fastq.gz | fastqc -q -t $ncores stdin -o fastqc_out/raw/combined
#rename the output from the previous command
mv fastqc_out/raw/combined/stdin_fastqc.html fastqc_out/raw/combined/combined_fastqc.html  
mv fastqc_out/raw/combined/stdin_fastqc.zip fastqc_out/raw/combined/combined_fastqc.zip  


#Assemble with PEAR
#This assumes your fastq files are named with [...]_1.fastq.gz or [...]_2.fastq.gz, where [...] is the sample ID
#Breakdown of the craziness below:
#  1) Find all gzipped FASTQ files and print their names
#  2) Use sed to only print out filename prefixes
#  3) Use sort to order the prefixes
#  4) Use uniq to remove duplicates (since fwd and rev files)
#  5) Take these files 4 at a time (since we have 4 cores)
#  6) For each file prefix, run PEAR in a subshell (multiprocessors)
#  7) Sleep for 60 seconds to let the 4 finish before starting next
#  8) Pipe the output to /dev/null so we don't get spammed
rm logPear.txt

numberOfFilePairs=$(ls $workingDir/sequence_files/*_1.fastq.gz| wc -l )
let numberOfIterations=($numberOfFilePairs+$ncores-1)/$ncores 

for j in $(seq 0 $numberOfIterations)
do
    let i=( $j * $ncores + 1 )
    echo $j "/" $numberOfIterations "..." 
    find $workingDir/sequence_files -name "*.fastq.gz" -printf '%f\n' | sed 's/_.*//' | sort | uniq | sed -n $i,$((i+${ncores}-1))p | while read line; do ( pear -f sequence_files/${line}_1.fastq.gz -r sequence_files/${line}_2.fastq.gz -o ${line} & ); done >> logPear.txt
    sleep 60
done

#Give a wee bit more time just in case any of the PEAR jobs have not caught up
#Sleep until we have
while [ $( ps -ef | grep pear | wc -l ) -gt 1 ]
do
 sleep 10
 echo "Waiting on PEAR to finish..."
done

# Get the percentage of reads assembled from the PEAR log files, check through these manually to make sure a high percentage of the pairs assembled. 
# If the percentages are low, there is a wide range of percentages, or there is one sample with an unusually low percentage, you should
# follow up by looking at FastQC results of the pre-merged files to figure out why this is the case before using any results.
grep "^Assembled reads" logPear.txt > logPearAssemblyStats.txt

#Ensure the sequence file is no longer here
rm -rf seq.fasta

#Take in the FASTQ files, funnel them into a single FASTA file
#FASTA headers in the form: ">SAMPLEID_SEQUENCENUMBER"
for filename in $( ls *.assembled.fastq )
do
    awk 'BEGIN{ORS=""; i=0;}{split(FILENAME, x, "."); prefix=x[1]; sub("@","",prefix); print ">" prefix "_" i "\n"; i+=1; getline; print; print "\n"; getline; getline;}' ${filename} >> seq.fasta
done

#Clean up our directory
mv *.fastq sequence_files
mkdir sequence_files/combined_fasta
mv seq.fasta sequence_files/combined_fasta/combined_seqs.fna


# Cluster using SortMeRNA and SUMACLUST via QIIME

inputFasta=$workingDir/sequence_files/combined_fasta/combined_seqs.fna

# prep param file
echo "pick_otus:threads " $ncores >> clustering_params.txt
echo "pick_otus:sortmerna_coverage 0.8" >> clustering_params.txt
echo "pick_otus:sortmerna_db $sortmernaDB" >> clustering_params.txt
echo "assign_taxonomy:id_to_taxonomy_fp $workingDir/reference_data/97_otu_taxonomy.txt" >> clustering_params.txt

# Use open reference clustering to first cluster by similarity to 16S reference sequences, 
# then use de novo clustering on the sequences that were not similar to anything in the database
# This will do clustering and taxonomic assignment
# full paths to all files must be provided
pick_open_reference_otus.py -i $inputFasta -o $workingDir/clustering/ -p $workingDir/clustering_params.txt -m sortmerna_sumaclust -s 0.1 -v --min_otu_size 1 --parallel --jobs_to_start $ncores  

# Remove low abundance OTUs 
# These may be due to cross-over from MiSeq lanes or sequencing errors OR they may be biological
# Make sure this is appropriate for your study before running
scripts/remove_low_confidence_otus.py -i clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o clustering/otu_table_high_conf.biom 

# Prepare summaries of the biom results to see the effect of removing low abundance OTUs
biom summarize-table -i clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o clustering/otu_table_mc1_w_tax_no_pynast_failures_summary.txt
biom summarize-table -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_summary.txt

#Convert the BIOM format file to tab-separated format (human readable)
biom convert --to-tsv -i clustering/otu_table_high_conf.biom --table-type='OTU table' -o clustering/otu_table_high_conf.tsv --header-key=taxonomy --output-metadata-id=Consensus\ Lineage

#We use awk on the converted OTU table to determine the lowest sequence depth
subsampleSize=$(awk 'BEGIN{FS="\t"} NR == 1 { } NR == 2 { max = NF-1; } NR > 2 { for (i = 2; i <= max; i++) { c[i] += $i; } } \
 END { smallest = c[2]; for (i = 3; i <= max; i++) { if (c[i] < smallest) { smallest = c[i]; }} print smallest; }' clustering/otu_table_high_conf.tsv)
echo $subsampleSize

#This is passed as a parameter to QIIME's rarefaction script
single_rarefaction.py -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_rarefied.biom -d $subsampleSize

#Convert rarefied table to tab-separated format
biom convert --to-tsv -i clustering/otu_table_high_conf_rarefied.biom --table-type='OTU table' -o clustering/otu_table_high_conf_rarefied.tsv --header-key=taxonomy --output-metadata-id=Cons$

### ANALYSIS OF OTU RESULTS ###

biomFile=clustering/otu_table_high_conf_rarefied.biom
biomTable=clustering/otu_table_high_conf_rarefied.tsv

#Taxonomy bar plots via QIIME
summarize_taxa_through_plots.py -f -s -i $biomFile -o taxaplot -m nwcs_arct_metadata.tsv

# Plot the beta diversity in 3D
beta_diversity_through_plots.py -m nwcs_arct_metadata.tsv -t clustering/rep_set.tre -i $biomFile -o qiime_pcoa_3D/

# To make 2D plots we can do the following
#Create distance matrices
beta_diversity.py -i $biomFile -o qiime_pcoa/distance/ -m weighted_unifrac -t clustering/rep_set.tre
beta_diversity.py -i $biomFile -o qiime_pcoa/distance/ -m bray_curtis

#Run PCoA on these matrices
principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa/

#Plot them
make_2d_plots.py -i qiime_pcoa/pcoa -m nwcs_arct_metadata.tsv -o qiime_pcoa/pcoa_plots/
chmod -R o+x qiime_pcoa/pcoa_plots


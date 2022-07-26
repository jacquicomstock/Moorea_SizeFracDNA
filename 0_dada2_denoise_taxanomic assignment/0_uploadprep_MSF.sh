# All the below commands were typed manually into the command line and not run as a script
# These are the commands I used to upload, move & unzip files, & run the other attached scripts as batch jobs in the cluster

#in local terminal
scp -r MiSeqAnalysis-20220720T224907Z-002.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Baratheon/
scp -r MiSeqAnalysis-20220720T224907Z-001.zip jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Baratheon/
scp -r silva_nr99_v138.1_wSpecies_train_set.fa.gz jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Baratheon/

#switch to cluster terminal window
#change to appropriate directory
cd Baratheon

#to unzip all .zip files
unzip '*.zip'

#move zipped files to a new folder to keep things organized
mkdir zipped_fastqs
mv *.zip zipped_fastqs/

#move SILVA database to the same directory as the unzipped fastqs
mv *.fa.gz MiSeqAnalysis/

#rename that folder to fastqs and move into that directory
mv MiSeqAnalysis fastqs
cd fastqs

#remove excess files imported (they were all zipped together)
rm -r Figures\ and\ working\ files/
rm -r filtered
rm -r *.csv
rm -r *.txt
rm -r *.jmp
rm -r *.xlsx
rm -r *.html
rm -r *.rds
rm -r *.xml
rm -r *.htm
rm -r *.Rmd

#make filtered directory for DADA2
mkdir filtered

#activate conda environment with R
conda activate R4.2.0

#make dada2 R script
vi Baratheon_dada2.R
#paste the dada2 R script into this new file and press esc
:w
:q

#make fasta R script
vi Baratheon_fasta.r
#paste the fasta R script into this new file and press esc
:w
:q

#submit batch job to slurm to run dada2 pipeline
sbatch \
	--job-name=Baratheon_dada2 \
	--nodes=1 \
	--tasks-per-node=32 \
	--cpus-per-task=1 \
	--mem=80G \
	--time=4:00:00 \
	--output=dada2_out \
	--error=dada2_err \
	--wrap="Rscript Baratheon_dada2.R"


#Once the R script has finished running, move output files to new folder
cd Baratheon
mkdir dada2_output_files
cd fastqs
mv *.fasta /home/jcomstock/Baratheon/dada2_output_files/
mv *.txt /home/jcomstock/Baratheon/dada2_output_files/
mv *.rds /home/jcomstock/Baratheon/dada2_output_files/
mv dada2_err /home/jcomstock/Baratheon/dada2_output_files/
mv dada2_out /home/jcomstock/Baratheon/dada2_output_files/

#Download the output files onto local computer, this is done in a local terminal window not in the cluster terminal window
scp -r jcomstock@pod.cnsi.ucsb.edu:/home/jcomstock/Baratheon/dada2_output_files/ /home/mobaxterm/Desktop/Research/Projects/Moorea/MOOREA\ 2019/Baratheon/0_CLEANED.REDONE.analyses/

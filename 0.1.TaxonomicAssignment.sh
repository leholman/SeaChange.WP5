#Here I assume you have blast (v2.12.0) & the nt database downloaded, mine is in the directory ~/data/gwm297.lholman/temp_db/nt and I will save the output to ~/data/gwm297.lholman/temp_db/nt/WP5.2/EUK|RIZ|MAM where the input fasta file is.  

#how to update the blast database with the perl script, these are already run so no need to do it again

# the script can be found in the below tar folder
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz 
#tar -xvf ncbi-blast-2.12.0+-x64-linux.tar.gz   

update_blastdb.pl --decompress nt 
update_blastdb.pl --decompress taxdb

#We first make a results folder
cd ~/data/temp_db/nt/WP5.2/RIZ|EUK|MAM
mkdir results 

#split fasta into chunks of 100 seqs for fast processing (do this for the 3 markers independently)
cd  ~/data/gwm297.lholman/temp_db/nt/
awk -v size=100 -v pre=split -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' file.fasta

##This is the raw expression we will use, the format listed is expected by other scripts/functions, dont mess with it.
#blastn -query WP5.1/split.00001 -db nt -out /WP5.1/results/test.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200

#Rather than submitting each expression we create a bunch of scripts that can then be submitted individually, there is probably a neater way to do this but I like this way as you can add many lines to be evaluated sequentially which I find hard to do with xsbatch

#make little scripts from within nt folder and give them execute permission
#submit scripts with 160000mb of RAM & 30 min limits
##RIZ

for i in split.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=160G\ncd ~/data/gwm297.lholman/temp_db/nt/\n./blastn -query WP5.2/RIZ/$bname -num_threads 1 -db nt -out WP5.2/RIZ/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > script.$bname.sh ;done
for i in ~/data/gwm297.lholman/temp_db/nt/WP5.2/RIZ/script.* ;  do  sbatch $i;done

##EUK
for i in split.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=160G\ncd ~/data/gwm297.lholman/temp_db/nt/\n./blastn -query WP5.2/EUK/$bname -num_threads 1 -db nt -out WP5.2/EUK/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > script.$bname.sh ;done
for i in ~/data/gwm297.lholman/temp_db/nt/WP5.2/EUK/script.* ;  do  sbatch $i;done

##MAM
for i in split.* ;do bname=$(basename $i); echo -e '#''!'"/bin/bash\n#SBATCH --ntasks=1\n#SBATCH --time=00:60:00\n#SBATCH --mem-per-cpu=160G\ncd ~/data/gwm297.lholman/temp_db/nt/\n./blastn -query WP5.2/MAM/$bname -num_threads 1 -db nt -out WP5.2/MAM/results/$bname.out -outfmt '6 qseqid qlen slen qcovs qcovhsp sseqid bitscore score evalue pident qstart qend sstart send staxids sscinames' -num_alignments 200" > script.$bname.sh ;done
for i in ~/data/gwm297.lholman/temp_db/nt/WP5.2/MAM/script.* ;  do  sbatch $i;done


#combine output from all the split files to form one mega taxonomy file
cd ~/data/temp_db/nt/WP5.2/results
cat split* > EUK.dada2.raw.taxonomy.txt
cat split* > RIZ.dada2.raw.taxonomy.txt
cat split* > MAM.dada2.raw.taxonomy.txt
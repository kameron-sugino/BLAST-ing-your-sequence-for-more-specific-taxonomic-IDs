-   [1 Introduction](#introduction)
    -   [1.1 Materials needed](#materials-needed)
        -   [1.1.1 BLAST database setup](#blast-database-setup)
        -   [1.1.2 Pulling required QIIME
            output](#pulling-required-qiime-output)
-   [2 Retrieving sequences of
    interest](#retrieving-sequences-of-interest)
-   [3 Running sequences through
    BLAST](#running-sequences-through-blast)
-   [4 Formatting the output](#formatting-the-output)

And then I started BLAST-ing: how to take unknown taxa sequences and
obtain a more specific ID

This is a short guide on using the NCBI BLAST command line software to
identify bacterial taxonomies using DNA sequences. For more info on
BLAST and installing BLAST+ see the following links (they’re not great
instructions, but will get you pointed in the right direction) \*
<https://www.ncbi.nlm.nih.gov/books/NBK52637/> \* Install BLAST+:
<https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/>

# 1 Introduction

-   The sequences used here are generated from a QIIME run of 16S V4-5
    sequences. The goal was to take unclassified taxa of interest (in
    this case, taxa unclassified to the order or class level for
    Clostridia/Clostridiales) and get a better ID on who is in the
    community

-   I ran the code on the Norman oscer servers, which had BLAST+
    downloaded for use. I only had to figure out how to download and
    unpack the reference database. The option “-remote” queues a run on
    the NCBI servers, so it’s not much different than BLAST-ing on their
    browser (i.e., fails if there are too many queries).

## 1.1 Materials needed

-   To start you will need
    -   BLAST database
    -   QIIME output tables

### 1.1.1 BLAST database setup

-   To download the database, follow the following:

-   Downloading 16S database on the server; all these steps occur on the
    shell, and should work if copy/pasted

    -   Make sure you know where all the files are being written! You’ll
        need to know their file location for later

``` r
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.??.tar.gz"

#extract contents of tarballs
for file in *.gz
do
tar -zxvpf "$file"
rm "$file"
done

wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"
tar -xzf taxdb.tar.gz
```

-   The untarred files are the database files you need. However, there
    are many “uncultured” and “environmental” sample names in the
    database, which are not particularly informative, so let’s remove
    them from consideration. Follow the steps below to reproduce the
    file.

-   I’d recommend remaking this file every so often just to make sure
    it’s the most recent version

-   Download and compile db for excluding any environmental, uncultured,
    unidentified, taxa

-   Use the subtree code from this github
    <https://github.com/pmenzel/taxonomy-tools#subtree>

    -   You’ll need to compile the code by downloading, upacking, and
        running “make” on the makefile
    -   You should be able to do this by navigating to the folder with
        the makefile:

``` r
cd /directory/with/makefile
make makefile
```

-   Download the taxa file to filter out undesired taxa assignments.
    Gotta FTP to the location (I used filezilla)
    <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>
    -   This file is provided in the data folder, but should be updated
        before each run to make sure it is the most recent version.
    -   I had to manually navigate to /pub/taxonomy/ and download the
        file, which is weird, but so it goes
-   Run the following to set up the exclusion list:

``` r
#unpack tarball
tar xf taxdump.tar.gz names.dmp nodes.dmp

#find names to exclude
grep "scientific name" names.dmp| grep -w -P "uncultured|environmental samples|metagenome|unclassified|unidentified" | cut -f1 >in.txt

#use subtree to make an exclusion list
/home/suginoka/Metagenomics_Programs/taxonomy-tools-master/src/subtree -t nodes.dmp -i in.txt > exclusion_list.txt
```

### 1.1.2 Pulling required QIIME output

-   Before running the program below on the HPC, you need to grab files
    from the QIIME output that end in the following:
    -   \*taxonomy.qza
    -   \*taxa-plot.qzv
    -   \*rep-seqs.qza
    -   \*table.qza
    -   \*table.qzv
-   These are basically zip files, so you can extract the contents into
    a folder using the native extract function in windows, with 7zip (or
    other compression/extraction program), or whatever the default
    option in macs is.
-   Within these unzipped folders, you’re looking for the following
    files (usually stored in the subfolder “data”):
    -   taxonomy.tsv
    -   level-1.csv through level-7.csv
        -   Don’t need these for this tutorial, but these are the
            kingdom-species OTU tables per sample
    -   dna-sequences.fasta
    -   feature-table.biom
    -   feature-frequency-detail.csv
-   The feature-table.biom file has to be processed into a tsv format,
    rather than a biom format.
    -   This table contains the OTU table corresponding to the seq_id
        (not taxa id, but the random string that tracks each unique
        sequence) per sample.
-   This conversion can be done with the following code and is
    relatively quick to run (note that the commands to convert to/from
    biom format are built into QIIME):

``` r
    biom convert -i feature-table.biom -o feature_otu_table_from_biom.txt --to-tsv
```

-   To find the sequences of interest, you’re going to need:
    -   taxonomy.tsv
    -   dna-sequences.fasta

# 2 Retrieving sequences of interest

-   Note that here I am specifically pulling Clostridia that are
    unclassified past the order or class levels
    -   The code will need to be rewritten to pull the specific taxa IDs
        of interest, so make sure you check that your output is correct

``` r
#merging taxonomy file with sequence file for BLAST analysis
tax<-read.table("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/taxonomy.tsv",sep="\t",header=T)
seq<-read.table("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/dna-sequences.fasta",header=F)

#reformat from fasta to table format
seq.df<-data.frame(Column1 = seq$V1[c(TRUE, FALSE)], 
                   Column2 = seq$V1[c(FALSE, TRUE)])
seq.df$Column1<-gsub(">","",seq.df$Column1)

#merge sequence table to taxa ID table
df<-merge(tax,seq.df,by.x="Feature.ID",by.y="Column1")

##################
#pull specific sequences of interest
df.1<-df[grep(c("Clostridia"),df$Taxon),]
#keep only unclassified order and class
#for many reasons, this isn't a great way of filtering taxa names, but I've basically capped the number of characters that can be in the string name since, if it's been classified past class, there will be more characters than 70 or so (hence choosing 75 as a cutoff)
df.o<-df.1[!nchar(df.1$Taxon)>75,]

#write.csv(df.o,"C:/Users/ksugino/Desktop/taxonomy_results/20230328_all_seqs_formatted.csv",row.names = F,quote=F)


#reformat into fasta format for BLAST
format<-data.frame()
a<-df.o
for(i in 1:nrow(a)){
  format<-rbind(format,paste0(">",a$Feature.ID[i]),a$Column2[i])
}

#write.table(format,"C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/202300801_clostridia_seqs.txt",row.names = F,quote=F)
```

# 3 Running sequences through BLAST

-   Now, you can put the files on the HPC and run something like the
    code below; note that you need your fasta-formatted sequences, the
    BLAST database, and the exclusion list to remove all
    uncultured/environmental IDs from the search.
    -   You’ll need to change the file locations for a few commands:
        -   Change /scratch/suginoka/BLAST_db/ to the location of the
            BLAST database
        -   Change 20230411_lacto_seqs.txt to the text file containing
            the fasta-formatted sequences
        -   Change /scratch/suginoka/BLAST_db/db/nt to the subfolder
            with the BLAST database (idk if this is needed, but it’s how
            I ran it and it works?)
        -   Change /scratch/suginoka/BLAST_db/exclusion_list.txt to the
            location/name of the file containing the list of taxa names
            to exclude from the search

``` r
    export BLASTDB=$BLASTDB:/scratch/suginoka/BLAST_db/

    module load BLAST+/2.13.0-gompi-2022a
    blastn -query 202300801_clostridia_seqs.txt -db /scratch/suginoka/BLAST_db/db/nt -out output.txt -outfmt "6 qseqid sseqid evalue pident stitle staxids sscinames scomnames sblastnames sskingdoms salltitles stitle" -evalue 1e-30 -task megablast -negative_taxidlist /scratch/suginoka/BLAST_db/exclusion_list.txt
```

# 4 Formatting the output

-   The next steps use the following files:
    -   output.txt (generated from the step above)
    -   feature-frequency-detail.csv
    -   The edited feature_otu_table_from_biom.txt output from above

``` r
df<-read.table("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/output.txt",fill=T,sep="\t",check.names = F,strip.white = F)

#this code is not the best and I'm so sorry (but it works!)

#if you see this an know of a better implementation (e.g., awk or something?) please let me know at kameron.sugino@gmail.com

#collect bacterial ID, and assembly score
query.id<-df[,1]
bac.name<-gsub("^([^ ]+ [^ ]+).*","\\1",df[,5])
bac.score<-df[,4]

#assemble and take the bacterial IDs with the max alignment score
#I don't know if there's a better way to do this, but here we are
reform<-data.frame(query.id,bac.name,bac.score)
reform.max<-aggregate(reform$bac.score, by = list(reform$query.id), max)
reform.max.merge<-merge(reform.max,reform,by.x=c("Group.1","x"),by.y=c("query.id","bac.score"))
reform.max.merge$bac.name<-gsub("\\[|\\]","",reform.max.merge$bac.name)

reform.table<-data.frame(table(reform.max.merge$Group.1,reform.max.merge$bac.name))
reform.table.e<-reform.table[reform.table$Freq>0,]

#read in and merge total counts (for all samples) for each seq ID
cnt<-read.csv("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/feature-frequency-detail.csv",header=F)

df.m<-merge(reform.table.e,cnt,by.x="Var1",by.y="V1")
colnames(df.m)<-c("Feature_ID","BLAST_Taxa_ID","BLAST_ID_Count","Sequence_Count")

final_ids<-aggregate(df.m$BLAST_ID_Count, by = list(df.m$Feature_ID), max)
final_ids.f<-merge(final_ids,df.m,by.x=c("Group.1","x"),by.y=c("Feature_ID","BLAST_ID_Count"))

#still have duplicate seq IDs, will combine names since it's unclear if there's a best bacterial ID
dupes<-final_ids.f[final_ids.f$Group.1 %in% final_ids.f$Group.1[duplicated(final_ids.f$Group.1)],]
dupes.m<-aggregate(dupes$BLAST_Taxa_ID,list(dupes$Group.1),paste, collapse = ",")

#create matrix of no dupes, of dupes, and then merge them
nodupes<-final_ids.f[!final_ids.f$Group.1 %in% final_ids.f$Group.1[duplicated(final_ids.f$Group.1)],]
colnames(nodupes)<-c("seq_id","BLAST_ID_Count","BLAST_Taxa_ID","Sequence_Count")

dupes.m.merge<-merge(dupes.m,dupes,by.x=c("Group.1"),by.y=c("Group.1"))
dupes.m.merge.e<-dupes.m.merge[!duplicated(dupes.m.merge$Group.1),]
colnames(dupes.m.merge.e)<-c("seq_id","BLAST_Taxa_ID","BLAST_ID_Count","Old_IDs","Sequence_Count")
dupes.f<-dupes.m.merge.e[,-4]

blast_final<-rbind(nodupes,dupes.f)

write.csv(blast_final,"C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/20230808_feature_summary_BLAST+_aligned.csv",row.names = F)
```

-   We still have some processing to do before the data are ready to
    use.

-   We’ll combine the feature summary generated above with the biom
    table and the metadata file (so we know which sample belongs to
    which treatment)

    -   Again, this code will likely need to be edited to fit your data

``` r
#take the new IDs and merge them with the feature otu table from the *.biom files
df<-read.csv("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/20230808_feature_summary_BLAST+_aligned.csv")
feat<-read.table("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/table.from_biom.txt",header=T,fill=T)
meta<-read.csv("C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/FMT-metadata.tsv",sep="\t")

#merge complete otu table with ids to calculate relative abundance later
otu.comp<-data.frame(gsub("X","",colnames(feat[,-1])),t(feat[,-1]))
otu.comp.f<-merge(meta,otu.comp,by.x="SeqID",by.y="gsub..X.......colnames.feat....1...")

final.otu<-merge(df,feat,by.x="seq_id",by.y="OTU")
#checking that the total number of reads pre seq_id is equivalent to the total number of reads imported from feature-frequency-detail.csv
rowSums(final.otu[,-c(1:4)])==final.otu$Sequence_Count

#need to merge taxa IDs again at this step, otherwise there will be multiple columns of taxa with the same bac. ID
otu.temp<-final.otu[,-c(1:2,4)]
otu.test<-aggregate(data = otu.temp, . ~ BLAST_Taxa_ID , sum)

id<-gsub("X","",colnames(final.otu[,-c(1:4)]))
otu<-data.frame(id,t(otu.test[,-c(1)]))
colnames(otu)<-c("id",otu.test$BLAST_Taxa_ID)

df.final<-merge(meta,otu,by.x="SeqID",by.y="id")

write.table(df.final,"C:/Users/ksugino/Desktop/Github_projects/BLAST-ing-your-sequence-for-more-specific-taxonomic-IDs/data/20230808_BLAST_OTU_table.txt",row.names = F, quote = F,sep='\t')
```

-   The output df.final is the otu file that contains the newly
    classified reads!
-   Note that this file only contains a *subset* of the total number of
    reads from microbiome, comprised of whatever sequences you pulled
    and BLASTED at the beginning

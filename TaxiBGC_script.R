#######################################################################################
## TaxiBGC (Taxonomy-guided Identification of Biosynthetic Gene Clusters): A computational pipeline for identifying experimentally verified Biosynthetic Gene Clusters (BGCs) and inferring their annotated SMs from metagenomic shotgun sequencing data. 
## The TaxiBGC pipeline includes three major steps: 
## 1) species-level, taxonomic profiling on the metagenome; 
## 2) a first-pass prediction of BGCs through querying the species (identified in the first step) in the TaxiBGC database; and 
## 3) confirmation (in silico) of the predicted BGCs (from the second step) based on the detection of BGC genes in the metagenome.

## Author: Utpal Bakshi, PhD and Vinod K. Gupta, PhD
## Date: September 2021

## Input: Paired-end fastq files (e.g., example_1.fq & example_2.fq)
## Output:'BGCs.csv' in current working directory

library(dplyr)
library(argparse)

#######################################################################################
## 1. Parsing the arguments from the user. Asks for paired-end fastq input files, directory of the TaxiBGC background database, and directory for the TaxiBGC output files
#######################################################################################
parser <- ArgumentParser()
parser$add_argument("-R1", "--R1Filename" , default=TRUE,
                    help="Path to input forward read file in fastq format")
parser$add_argument("-R2", "--R2Filename" , default=TRUE,
                    help="Path to input reverse read file in fastq format")
parser$add_argument("-o", "--directory" , default=TRUE,
                    help="Name of output directory to be created")
parser$add_argument("-d", "--database" , default=TRUE,
                    help="Path to TaxiBGC databbase directory")
args <- parser$parse_args()
input1 <- args$R1Filename
input2 <- args$R2Filename
database <- args$database
sample_name <- args$directory
dir.create(sample_name)
path1<-paste(" ",sample_name,'/sampleID_1',sep="")
path2<-paste(" ",sample_name,'/sampleID_2',sep="")
print(path1)
system(paste('cp ', input1,path1))
system(paste('cp ', input2,path2))

#######################################################################################
## 2. Running MetaPhlAn2 on metagenome
#######################################################################################
setwd(sample_name)
dir.create("temp")
system('metaphlan sampleID_1,sampleID_2 --bowtie2out temp/sampleID.bowtie2.bz2 --input_type fastq > temp/MetaPhlAn2_raw.txt')

#######################################################################################
## 3. Generating species abundance matrix from MetaPhlAn2 output
#######################################################################################
setwd("./temp")
system('awk -F"\t" \'/__|clade_name/ {print $1"\t"$3}\' example.tsv > MetaPhlAn2.txt')
system('grep -E \"(s__)|(^ID)\" MetaPhlAn2.txt | grep -v \"t__\" | sed \'s/^.*s__//g\' > ./Species_MetaPhlAn2.txt')
system('awk \'NR==1 { print \"#SampleID\", FILENAME }; 1\' \"Species_MetaPhlAn2.txt\" > temp && mv temp \"Species_MetaPhlAn2_with_header.txt\"') # Adding filenames in file header

#######################################################################################
## 4. Analysis of MetaPhlAn2 output
#######################################################################################
df <- read.csv('Species_MetaPhlAn2_with_header.txt', header= TRUE, sep = "\t",  check.names = F, row.names = 1)
df1 <- cbind(rownames(df), data.frame(df, row.names=NULL))
remove <- c('noname_', 'unclassified','virus') ##  remove 
modified_df1 <- df1[ !grepl(paste(remove, collapse="|"), df1$`rownames(df)`),]
modified_df1 <- data.frame(modified_df1[,-1], row.names = modified_df1[,1])
df2<-sweep(modified_df1,2,colSums(modified_df1),'/') # re-normalize, sum of a sample's relative abundances becomes 1
df2[is.na(df2)] <- 0
df2[df2 < 0.00001] <- 0 # criteria for species being present is >= 0.001% 
df3 <- cbind(rownames(df2), data.frame(df2, row.names=NULL))
colnames(df3) <- c("Species_name","Abundance")
df4 <- df3[df3$Abundance != 0, ]

#######################################################################################
## 5. Overlap with background database Step1
#######################################################################################
genome_file <- read.csv(file.path(database,"TaxiBGC_database.csv"), header= TRUE, sep = ",", check.names = F)
df2_match<-merge(as.data.frame(df4),as.data.frame(genome_file), by.x = 'Species_name',by.y = 'Species_name', all=F)
df2_bT2<-df2_match[, colSums(df2_match != 0) > 0] # remove all blank samples
## To get step 1 results
df2_bT2$products <- paste(df2_bT2$BGC_Accession,df2_bT2$Species,df2_bT2$Product)
products_subset<-subset(df2_bT2, !duplicated(subset(df2_bT2, select=c(products))))
products_final<-products_subset[, colSums(products_subset != 0) > 0] # remove all blank samples
Step1_results<- data.frame(subset(products_final, select = -c(products)))
Step1_results <- subset(Step1_results, select = -c(Abundance, Similarity))

#######################################################################################
## 6. Fetch BGCs
#######################################################################################
bgc_name <- data.frame(unique(Step1_results$BGC_Accession))
write.table(bgc_name, file = "bgc_name.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)

#######################################################################################
## 7. Preparing files for bowtie2
#######################################################################################
dir.create ("./temp2")
system('find . -type f ! -name \'bgc_name.txt\' -delete')
### Copy bgc fasta 
BGC_fasta_path<-file.path(database,"bgc_fasta","*")

system(paste('cp ',BGC_fasta_path,' .'))
system('for file in *.txt; do xargs < $file cat > ./temp2/$file; done')
file.copy("./temp2/bgc_name.txt", "./../")
setwd("./..")

#######################################################################################
## 8. Running bowtie2 
#######################################################################################
system('bowtie2-build bgc_name.txt bgc_seq')
system('bowtie2 -x bgc_seq -1 sampleID_1 -2 sampleID_2 > bgc_seq.bam')
system('samtools sort -o sorted.bam bgc_seq.bam')
system('samtools index sorted.bam')
system('samtools idxstats sorted.bam > ./temp/bgc_idxstats.txt')

#######################################################################################
## 9. Cleaning bowtie2 results
#######################################################################################
setwd("./temp")
system('awk \'{ print $1, $2, $3 + $4; }\' bgc_idxstats.txt > bgc_idxstats.fasta') # Idxstats to count
system('awk \'BEGIN{print \"contig length bgc_idxstats.txt\"}1\' bgc_idxstats.fasta > bgc_idxstats2.fasta') # Adding headers
system('cat bgc_idxstats2.fasta | tr \' \' \'\t\' > bgc_idxstats3.fasta') # Tab delimitted
system('awk \'{sub(/_.*/, \"\", $1)} 1\' bgc_idxstats3.fasta > bgc_idxstats4.txt') # Trimming first column
system('awk \'{print $1 \" \" $3}\' bgc_idxstats4.txt > bgc_idxstats5.txt') # Removing second column
system('find . -type f ! -name \'bgc_idxstats5.txt\' -delete')

#######################################################################################
## 10. Processing of bowtie2 results
#######################################################################################
filenames <- list.files(full.names=F, pattern=".txt")
output <-lapply(filenames,function(i){
  t<-read.csv(i, header=T, check.names = F, sep = " ")
  if(all(t[[2]] == 0)) return(NULL)
  t$gene_count<-1
  t[,2][t[,2]>0]<-1
  presence_absence_df<-aggregate(. ~ contig, t, sum)
  presence_absence_df$sample_name<-names(t[2])
  colnames(presence_absence_df)<-c("BGC_Accession","Gene_presence", "Gene_count", "Sample_name")
  presence_absence_df$Percentage<-(presence_absence_df$Gene_presence/presence_absence_df$Gene_count)*100
  presence_absence_df<-presence_absence_df[presence_absence_df$Percentage != 0, ]
  presence_absence_df$tp_step2_20_percent<-length(presence_absence_df$Percentage[presence_absence_df$Percentage>=20])
  presence_absence_df<-presence_absence_df[presence_absence_df$Percentage >= 20, ]
  presence_absence_df <- subset(presence_absence_df, select = -c(Gene_presence, Gene_count, Percentage))
  colnames(presence_absence_df)<-c("BGC_name", "Sample", "BGCs_step2_20_percent")
  presence_absence_df <- presence_absence_df [c("Sample", "BGCs_step2_20_percent", "BGC_name")]
})
Step2_results2_20<-do.call(rbind,output)
Step2_results2_20<-subset(Step2_results2_20, BGC_name!="*")

#######################################################################################
## 11. Rearrange the output files
#######################################################################################
Step2_results2_20[,1] <- data.frame(gsub("_idxstats.txt.*$", "", Step2_results2_20[,1]))
Step2_results_20_matrix <-tidyr::pivot_wider(Step2_results2_20, names_from = Sample, values_from = BGCs_step2_20_percent)
Step2_results_20_matrix[is.na(Step2_results_20_matrix)] <- 0
Step2_results_20_matrix[-1] <- + sapply(Step2_results_20_matrix[-1], as.logical)
Step2_results_20_matrix<-Step2_results_20_matrix[apply(Step2_results_20_matrix[,-1], 1, function(x) !all(x==0)),] # Removing all zero rows
Step2_results_20_matrix<-Step2_results_20_matrix[, colSums(Step2_results_20_matrix != 0) > 0] # Removing all zero columns
prevalence_20_matrix<- (rowSums(Step2_results_20_matrix[,-c(1)]!=0)/length(Step2_results_20_matrix[,-c(1)]))*100
Step2_results_20_matrix<-cbind(Step2_results_20_matrix,prevalence_20_matrix)
Step2_results_20_matrix<- subset(Step2_results_20_matrix, select = -c(bgc,prevalence_20_matrix))
Step2_results<-merge(as.data.frame(Step2_results_20_matrix),as.data.frame(Step1_results), by.x = 'BGC_name',by.y = 'BGC_Accession', all=F)
colnames(Step2_results)[1] <- "BGC_Accession"

#######################################################################################
## 12. Final clean-up
#######################################################################################
setwd("./..")
unlink("temp", recursive = TRUE)
unlink("*.bt2", recursive = TRUE)
unlink("*.bam", recursive = TRUE)
unlink("*.sh", recursive = TRUE)
unlink("*.txt", recursive = TRUE)
unlink("*.bai", recursive = TRUE)
unlink("sampleID_1", recursive = TRUE)
unlink("sampleID_2", recursive = TRUE)
write.csv(Step2_results,"BGCs.csv")


library(gwascat)
#for large text files, use fread in library data.table (or the older version read.table) is more efficient
#or second option is using read_table in readr, slower than fread in data.table
#third option: read.csv.raw from iotools for quickly reading csv files.

library(data.table)
#setwd("C:/Users/yuy/Desktop/kaggle/DNA") , set your own path
mydata<-fread('MyDNA.txt')

library(ggplot2)
mydata$chromosome = ordered(mydata$chromosome, levels=c(seq(1, 22), "X", "Y", "MT"))
ggplot(mydata, aes(chromosome)) + geom_bar()+labs(title="Bar plot of snps on chromosome" )

#use hg18 because 23andMe used this transcript annotation database as reference
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
objects(txdb)
#txdb offers GenomicRange objects! 

#group annotated transcripts  by  genes, other options are grouping by CDS and exons (cdsBy, exonsBy)
tx.by.gene <- transcriptsBy(txdb,"gene")

#convert the Entrez Gene IDs from txdb to actual gene names, gene name symbol is the "key" field
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
select(org.Hs.eg.db, keys="ADH1B", columns=c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL")

#build a GRanges object from my genotyping data, and look for overlaps
levels(mydata$chromosome) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")

my.snps <- with(mydata, GRanges(seqnames=chromosome, 
                           IRanges(start=position, width=1), 
                           rsid=`# rsid`, genotype=genotype)) 

#determine overlaps between my genotype data and the tx.by.gene
#the IRanges locations are site of mutations.
adh.i <- findOverlaps(tx.by.gene["125"], my.snps)
hits <- subjectHits(adh.i)
my.snps[hits]
#note, Herve Pages changed the findOverlaps value. It is now an object of the Hits class that does not support the matchMatrix accessor anymore. The code was adapted to the new accessor subjectHits.


##################################### look for mutations using gwascat ##########################

#first load the metadata gwrngs19 for human DNA
gwrngs.emd <- as.data.frame(elementMetadata(gwrngs19))
dm <- merge(mydata, gwrngs.emd, by.x="# rsid", by.y="SNPs")
objects(dm)

################################### Look at high risk alleles #############################
risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)
#use mapply to vectorize the output with multivariate input
#the code below takes vector risk.allels and check for them in dm$genotype which is my genotype.
my.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))
}, risk.alleles, dm$genotype)

dm$my.risk <- my.risk
risk_cal<- dm[dm$my.risk, ]

rel.cols <- c(colnames(mydata), "Disease.Trait", "Risk.Allele.Frequency",
              "p.Value", "my.risk", "X95..CI..text.")
ordered_data<-order(risk_cal$Risk.Allele.Frequency)
select_data<-subset(risk_cal[ordered_data], select=rel.cols)
head(select_data)[1]
#find index corresponding to hypothyroidism in disease trait
idx_thyroid<-grep("Hypothyroidism", select_data$Disease.Trait)
thyroid_rsid<-select_data[idx_thyroid]$`# rsid`
#match the rsid to dm  
thyroid_data<-dm[match(thyroid_rsid, dm$`# rsid`)]
head(thyroid_data)
#find studies that are specific to Chinese
idx<-grep("chinese", thyroid_data$Initial.Sample.Size)
head(select_data[idx])


############################### Plot the karyograms #################################
library(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
head(hg19IdeogramCyto)

w <- plotStackedOverview(hg19IdeogramCyto, cytoband=FALSE)
(elementMetadata(gwrngs19)$my.genotype <- 
    mydata$genotype[(match(elementMetadata(gwrngs19)$SNPs, mydata$`# rsid`))])

elementMetadata(gwrngs19)$risk_cal<- with(elementMetadata(gwrngs19), 
                                        mapply(function(risk, mine) {
                                          risk %in% unlist(strsplit(mine, ""))
                                        }, gsub("[^\\-]*-([ATCG?])", "\\1", Strongest.SNP.Risk.Allele), my.genotype))

p + layout_karyogram(gwrngs19, aes(color=risk_cal))

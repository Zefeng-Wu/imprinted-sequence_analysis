library(rtracklayer)
library(GenomicFeatures)
library(Biostrings)
data<-read.table("1out",stringsAsFactors = FALSE)
repeat_gr<-makeGRangesFromDataFrame(data,keep.extra.columns = TRUE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")

repeat_gr$V15<-DNAStringSet(repeat_gr$V15) # take seed sequence as biostring
repeat_gr$V16<-DNAStringSet(repeat_gr$V16) #take seed sequence as biostring
repeat_gr$repeat_gc<-rowSums(letterFrequency(repeat_gr$V15,"CG",OR=0,as.prob = TRUE))


tr<-makeTxDbFromGFF("~/MyResearch/genome.db/TAIR/gtf/Arabidopsis_thaliana.TAIR10.31.gtf")
genes<-genes(tr)
imp_gene<-read.table("~/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/ara_imp_by_paper_num/imp2+.list",stringsAsFactors = FALSE) # high confidence
imp_gr<-genes[names(genes)%in%imp_gene$V1]

## repeat_gr overlapped with imp (gene body)
length(unique(queryHits(findOverlaps(query = imp_gr,subject = repeat_gr))))
#length(unique(queryHits(findOverlaps(query = imp_gr,subject = repeat_gr)))) #13 imprinted genes contain repeat
length(unique(queryHits(findOverlaps(query = genes,subject = repeat_gr))))   #2658
phyper(13-1,2658,33602-2658,90,lower.tail = FALSE) # 0.024
overlaped_gene_repeat<-repeat_gr[subjectHits(findOverlaps(query = imp_gr,subject = repeat_gr))] # overlaped repeat grange

## length distribution 
library(ggplot2)
df1<-data.frame(length=width(overlaped_gene_repeat$V16),class=rep("Imp",length(overlaped_gene_repeat)))
df2<-data.frame(length=width(repeat_gr$V16),class=rep("All",length(repeat_gr)))
df<-rbind(df1,df2)
ggplot(df,aes(x=class,y=log(length)))+geom_boxplot()+ylab("Length(log10)")

## gc content
df1<-data.frame(gc= overlaped_gene_repeat$repeat_gc,class=rep("Imprinted ",length(overlaped_gene_repeat)))
df2<-data.frame(gc= repeat_gr$repeat_gc,class=rep("All",length(repeat_gr)))
df<-rbind(df1,df2)
ggplot(df,aes(x=class,y=gc))+geom_boxplot()+ylab("GC content")+theme(text = element_text(size = 25))


## repeat_gr overlapped with imp (up1000 with gene body)
gene_with_up1000_gr<-punion(flank(imp_gr,width = 1000),imp_gr)
findOverlaps(query = gene_with_up1000_gr,subject = repeat_gr) #18 imprinted genes

all_genes_with_up_1000 <- punion(flank(genes,width = 1000),genes)  
findOverlaps(query = all_genes_with_up_1000,subject = repeat_gr)  #4490 genes
phyper(18-1,4490,33602-4490,90,lower.tail = FALSE) # 0.05

#############
## 525 ara imp
#'''
imp_gene<-read.table("~/MyResearch/Imprinting_prediction/imprint_gene_list/TAIR10/6imprinted.list",stringsAsFactors = FALSE)  # low confidence
imp_gr<-genes[names(genes)%in%unique(imp_gene$V1)]

## overlap
length(unique(queryHits(findOverlaps(query = imp_gr,subject = repeat_gr))))  #55 imprinted genes contain repeat
length(unique(queryHits(findOverlaps(query = genes,subject = repeat_gr))))   #2658
phyper(55-1,2658,33602-2658,525,lower.tail = FALSE) # 0.02

## neary by 1000
gene_with_up1000_gr<-punion(flank(imp_gr,width = 1000),imp_gr)
findOverlaps(query = gene_with_up1000_gr,subject = repeat_gr) #92 imprinted genes

all_genes_with_up_1000 <- punion(flank(genes,width = 1000),genes)  
findOverlaps(query = all_genes_with_up_1000,subject = repeat_gr)  #4490 genes
phyper(92-1,4490,33602-4490,525,lower.tail = FALSE) # 0.003


#### repeat distribution (soggi)
library(Biostrings)
fa<-readDNAStringSet("~/MyResearch/genome.db/TAIR/dna/Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa")

repeat_gr<-makeGRangesFromDataFrame(data,keep.extra.columns = FALSE,ignore.strand = TRUE,seqnames.field = "V1",start.field = "V2",end.field = "V3")
seqlengths(repeat_gr)<-width(fa)[match(names(seqlengths(repeat_gr)),names(fa))] 

reduce_gr<-reduce(repeat_gr)
reduce_gr$score<-1
export.bw(reduce_gr,"out.bw") # exist overlap between repeat, so reduce

library(soGGi)
hipExample1 <- regionPlot("out.bw",testRanges =genes,format = "bigwig",distanceUp = 1000,distanceDown = 1000,style = "percentOfRegion",method = "spline")
hipExample2 <- regionPlot("out.bw",testRanges =imp_gr,format = "bigwig",distanceUp = 1000,distanceDown = 1000,style = "percentOfRegion",method = "spline")

gl<-GRangesList(genes,genes[genes$gene_id%in%imp_gr$gene_id])
names(gl)<-c("All","Imp")
plotRegion(hipExample1,gts = gl)  # arrary 1*2 plot

#### self plot for soggi
require(TTR)
data_a = plotRegion(hipExample1)$data
data_b = plotRegion(hipExample2)$data
data_all<-data.frame(All=data_a$Score,Imp=data_b$Score,Region= seq(1:300))

data_all$Imprinted <- SMA(data_all$Imp, 20)
ggplot(data_all,aes(Region))+geom_line(aes(y = All, colour = "All"))+
                            geom_line(aes(y =Imprinted, colour = "Imprinted"))+
                            scale_x_discrete(limits=seq(1,400,100),labels=c("1000","TSS","TTS","1000"))+
                            ylab("Tandem repeat density")+
                            labs(colour = "Class")
  

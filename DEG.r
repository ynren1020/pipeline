#####################2019-07-23####################
##differential expression analysis#################
##then GSEA preranked analysis#####################
##pipeline from DEA and rank file prep#############
###################################################

####load libraries####
library("edgeR")
#library("sva")
#library("R.utils")
library("limma")


#input<-"DMS53_shASCL1.txt"
#input<-"HUVEC.txt"
#input<-"DMS_53_Luc_DARPP_OE.txt"
DEG("HUVEC.txt",6,9,2,2) #test if this function work#
DEG("DMS_53_Luc_DARPP_OE.txt",6,11,3,2)

DEG<-function(input,col_start,col_end,trials,times){

dat <- read.delim(input, row.names="Geneid", check.names=FALSE,skip=1);
datsub <-dat[,col_start:col_end]
trials<-trials
group <- factor(rep(1:trials,each=times))


y <- DGEList(counts=datsub, group=group);
y$genes <- data.frame(Length=dat$Length);
rownames(y$genes) <- rownames(y$counts);

#y$genes$Symbol <- data.frame(Length=src.data$Symbol);
#rownames(y$symbol) <- rownames(y$counts);

# filtering out low expressed genes
#  the minimum number of samples in each group is 2, over here.
keep <- rowSums(cpm(y)>1) >= times
table(keep);
# FALSE  TRUE 
# 41526 15607
y <- y[keep, keep.lib.sizes=FALSE];



design <- model.matrix(~0+group);
rownames(design) <- colnames(y)
design


# Note that the filtering does not use knowledge of what treatment corresponds to each sample, so
# the filtering does not bias the subsequent differential expression analysis.
# The TMM normalization is applied to account for the compositional biases:


# normalization by the library sizes
y <- calcNormFactors(y);
y$samples;

write.table(rpkm(y), file=paste0(strsplit(input,"[.]")[[1]][1],".rpkm.count.txt"), sep="\t", row.names = TRUE, quote = FALSE);


#barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
#title("Barplot of library sizes")

# normalise the read counts with "voom" function
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

#boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let"s add a blue horizontal line that corresponds to the median logCPM
#abline(h=median(v$E),col="blue")

# save normalised expression data into output dir
write.table(counts.voom,file=paste0(strsplit(input,"[.]")[[1]][1],".voom.count.txt"),row.names=T,quote=F,sep="\t");


########################################################################################################################
# fit linear model for each gene given a series of libraries
for (i in 2:trials)
  {
fit <- lmFit(v, design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
#matrix.i<-paste0("matrix",i,"vs1")

matrix.i <- makeContrasts(contrasts = paste0("group",i,"-group1"),levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
#fit.i<-paste0("fit",i,"vs1")
fit.i<- contrasts.fit(fit, matrix.i)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.i <- eBayes(fit.i)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.i, p.value=0.05,lfc=1))
# group2 - group1
# Down               0        
# NotSig           13236      
# Up                 0        

num = length(fit.i$genes$Length)
#degs.i<-paste0("degs",i,"vs1")
degs.i <- topTable(fit.i, coef=paste0("group",i,"-group1"), confint=TRUE, number = num)
write.table(degs.i, file=paste0(strsplit(input,"[.]")[[1]][1],".degs.test.txt"), sep="\t",row.names = TRUE, quote = FALSE);

##make preranked list for GSEA##
#rank.i<-paste0("rank",i,"vs1")
rank.i<-data.frame("gene_name"=row.names(degs.i),"rank"=degs.i$logFC)
write.table(rank.i, file=paste0(strsplit(input,"[.]")[[1]][1],paste0("rank",i),".vs1.test.txt"), sep="\t",row.names = FALSE, col.names=FALSE,quote = FALSE);

}

}






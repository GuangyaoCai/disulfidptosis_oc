setwd("D:/工作妇肿/课题/卵巢_预后模型v2/TCGA_ov/")
exp <- read.table(file ="TCGA-OV.htseq_counts.tsv.gz", header=T,check.names=F, sep = "\t") #exp格式是log2(count+1)
head(exp[1:4,1:4])

names(exp)[1]<-"id" 
exp$id[which(exp$id == "AES")] <- "TLE5"
rownames(exp)<- exp$id
head(exp[1:4,1:4])

exp[1:4, 1:4]
library(edgeR)
range(exp[,-1])
cex<-cpm(2^(exp[,-1])-1)#转化为cpm格式
cex<-log2(cex+0.1)
range(cex)
cex1<- cbind(id=rownames(cex),cex)
cex1 <- t(cex1)
cex1[1:4, 1:4]
dim(cex1)

cex2 <- cex1[-1,]
cex2 <- cbind(sample= rownames(cex2),cex2)
cex2[,1]<-substr(cex2[,1],1,12)
cex2[1:4, 1:4]
dim(cex2)

library(limma)
library(magrittr)
setwd("D:/工作妇肿/课题/卵巢_预后模型v2/TCGA_ov/")
cex3<- cex2

row.names(cex3)<- substr(rownames(cex3),1,12)
dim(cex3)
table(duplicated(rownames(cex3)))
cex3[1:4, 1:4]
range(cex3[,-1])
cex3<-avereps(as.matrix(cex3[,-1]))
head(cex3[1:4,1:4])
dim(cex3)
cex4 <- t(cex3) %>% cbind(id= row.names(.), .)
dim(cex4)
cex4[1:4, 1:4]

write.table(cex4,file = "log2cpm___new.txt", sep="\t",row.names=F,quote=F)

rm(list=ls())

setwd("D:/工作妇肿/课题/disulfidptosis_v2/")
gene<-read.table("disulfidptosis.txt",check.names = F,header = F)
colnames(gene)<- c("gene")

setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
rt<-read.table("log2cpm___new.txt",check.names = F)
head(rt[1:4,1:4])
row.names(rt) <- rt[,1]
colnames(rt) <- rt[1,]
rt1 <- t(rt)[-1,]
head(rt1[1:4,1:4])
colnames(rt1)[1]<- c("id")
geneList <- colnames(rt1)[-1]

geneList <- sapply(strsplit(geneList, "\\."), "[", 1)

library(clusterProfiler)

gene.df <- bitr(geneList, fromType="ENSEMBL",
                toType="SYMBOL", 
                OrgDb = "org.Hs.eg.db")
gene.df <- na.omit(gene.df)

rt2 <- t(rt1) %>% as.data.frame() %>% cbind(ENSEMBL=rownames(.), .)

rt2[1:4, 1:4]
rt2$ENSEMBL <- sapply(strsplit(rt2$ENSEMBL, "\\."), "[", 1)

rt3<- merge(gene.df, rt2, by= "ENSEMBL")
rt3[1:4, 1:4]
rt3<- as.matrix(rt3)
row.names(rt3)<- rt3[,2]
rt3<- rt3[,-1]
rt3[1:4, 1:4]

library(limma)
library(magrittr)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")

rt4<-avereps(as.matrix(rt3[,-1]))
head(rt4[1:4,1:4])
dim(rt4)
total_gene <- rownames(rt4)


exp1<-rt4[which(rownames(rt4)%in%gene$gene),]  %>% cbind(id= rownames(.), .) %>% as.data.frame()
setdiff(gene$gene, exp1$id)

survival<-read.table("all_survival.txt",header=T,sep="\t",check.names=F,row.names=1)
head(survival[1:4,1:2])
survival1 <- cbind(rownames(survival),survival)
survival1[,1]<-substr(survival1[,1],1,12)
colnames(survival1)[1]<- c("id")
survival1 <- unique(survival1[,1:3])
rownames(survival1) <- survival1$id

setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
exp2 <- t(exp1[,-1])
exp2 <- cbind(id=row.names(exp2), exp2)
final_sur_data<-merge(survival1, exp2, by="id")
head(final_sur_data[1:4,1:4])

rownames(final_sur_data)<- final_sur_data$id
head(final_sur_data[1:4,1:4])
write.table(final_sur_data,file = "disulfidptosis_sur_exp_new.txt",sep="\t",row.names=F,quote=F)


rm(list = ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
final_sur_data <- read.table("disulfidptosis_sur_exp_new.txt", sep="\t", header = T, row.names = 1)
final_sur_data_OS1 <- final_sur_data
head(final_sur_data_OS1[1:4,1:4])
sampleID <- row.names(final_sur_data_OS1)
final_sur_data_OS1 <- sapply(final_sur_data_OS1, as.numeric)
final_sur_data_OS1 <- as.data.frame(final_sur_data_OS1)
head(final_sur_data_OS1[1:4,1:4])
row.names(final_sur_data_OS1) <- sampleID
head(final_sur_data_OS1[1:4,1:4])
final_sur_data_OS1$id <- rownames(final_sur_data_OS1)

write.table(final_sur_data_OS1, "final_sur_data_OS1_new.txt", sep="\t",row.names=F,quote=F)

pFilter=1
outTab=data.frame()
sigGenes=c("fustime","fustat")
for(gene in colnames(rt[5:ncol(rt)])){
  diff=survdiff(Surv(fustime, fustat) ~rt[[gene]],data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  #单因素cox分析
  cox=coxph(Surv(fustime, fustat) ~ rt[[gene]], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(gene=gene,
                     KM=pValue,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxP) )
  if((pValue<pFilter) & (coxP<pFilter)){
    sigGenes=c(sigGenes,gene)
  }
}
setwd("D:/工作妇肿/课题/disulfidptosis_v2/OS_p0.05/")
write.table(outTab,file = "2_cluster_survival_HR_surv_OS_new.txt",row.names = T)


#3.LASSO回归精简变量

library("glmnet")
library("survival")
setwd("D:/工作妇肿/课题/disulfidptosis_v2/OS_p0.05/")
la<- read.table("2_cluster_survival_HR_surv_OS_new.txt")

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
rownames(final_sur_data_OS1) <- final_sur_data_OS1$id
final_sur_data_OS2 <- t(final_sur_data_OS1)
head(final_sur_data_OS2[1:4,1:4])
final_sur_data_OS3 <- cbind(gene=rownames(final_sur_data_OS2), final_sur_data_OS2)
la1 <- merge(la, final_sur_data_OS3, by = "gene")
la3 <- t(la1)
pvalue1<- la3[2,]
colnames(la3)<- la3[1,]
la4<- cbind(final_sur_data_OS1[,c(1,2)],la3[-c(1:6),])
head(la4[1:4,1:4])
la4<-avereps(as.matrix(la4))
head(la4[1:4,1:4])
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
write.table(la4,file="LASSO_input.txt",sep="\t",quote=F,row.names=F)

set.seed(123456)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
rt=read.table("LASSO_input.txt",header=T,sep="\t")     
head(rt[1:4,1:4])

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$fustime,rt$fustat))

fit <- glmnet(x, y, family = "cox", alpha=1, maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene1 <- outTab[which(outTab$gene %in% lassoGene),]
write.table(lassoGene1,file="lassoGene.txt",sep="\t",quote=F,row.names=F,col.names=T)


setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
gene<-read.table("lassoGene.txt",header=F,sep="\t",check.names=F)
rownames(gene)<-gene$V1
colnames(gene)<- gene[1,]
sur<-final_sur_data_OS1[,c(1,2)]
exp<-final_sur_data_OS1[,which(colnames(final_sur_data_OS1)%in%rownames(gene))]

mc<-cbind(sur,exp)
mc<-cbind(id=rownames(mc),mc)
head(mc[1:4,1:4])
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
write.table(mc,file="laOutputExp.txt",sep="\t",quote=F,row.names=F)

#4.multicox
library(survival)
library(survminer)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
gene<-read.table("lassoGene.txt",header=F,sep="\t",check.names=F)
rownames(gene)<-gene$V1
dim(gene)

sur<-final_sur_data_OS1[,c(1,2)]
exp<-final_sur_data_OS1[,which(colnames(final_sur_data_OS1)%in%rownames(gene))]

mc<-cbind(sur,exp)
mc<-cbind(id=rownames(mc),mc)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
write.table(mc,file="multiInput.txt",sep="\t",quote=F,row.names=F)

####multicox
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
multiInput=read.table("multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
head(multiInput[1:4,1:4])
dim(multiInput)

sampleID <- rownames(multiInput)

multiCox=coxph(Surv(fustime, fustat) ~ ., data = multiInput)
multiCox=step(multiCox,direction = "forward")
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.csv(outTab,file="multiCox.csv",row.names=F,quote=F)


library(survival)
library(survminer)
library(eoffice)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_lasso/")
pdf(file = "multiCox_forest.pdf",width=7,height=4, onefile = F)
ggforest(model=multiCox,
         data = multiInput,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.5, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

riskScore=predict(multiCox,type="risk",newdata=multiInput)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("fustime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"High_risk","Low_risk"))
library(magrittr)
multiInput_2 <- cbind(multiInput[,outCol],riskScore,risk) %>% cbind(id=rownames(.), .)
median(riskScore) ####riskscore的cutoff值

risk <- multiInput_2 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
write.table(multiInput_2,file="risk_train.txt", sep="\t",quote=F,row.names=F)

#riskScore
library(survival)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
multiInput_2=read.table("risk_train.txt",header=T,sep="\t")
diff=survdiff(Surv(fustime,fustat) ~risk,data = multiInput_2)
pValue=1-pchisq(diff$chisq,df=1)
#pValue=round(pValue,3)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(fustime,fustat) ~ risk, data = multiInput_2)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
pdf(file="survival_train.pdf",width=5.5,height=5, onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata", 
           linetype = "strata",
           surv.median.line = "hv", 
           ggtheme = theme_bw(), 
           palette = c( "#E4392E","#2546F1"),
           xlim= c(0, 156),
           break.x.by= 24)
dev.off()

multiInput_2$risk <- as.factor(multiInput_2$risk)
a <- c("Low_risk", "High_risk")
multiInput_2$risk <- factor(multiInput_2$risk, levels=a)
riskHR <- coxph(Surv(multiInput_2$fustime,multiInput_2$fustat)~risk, data=multiInput_2)
riskHR

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
save(multiCox, survival, coxGene, file = "multiCox.RData")


#验证
rm(list = ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
load("multiCox.RData")

outCol=c("fustime","fustat",coxGene)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GEO/validation_data/")

library(magrittr)
library(survival)
library(survminer)

#E_MTAB_386
load("D:/工作妇肿/课题/disulfidptosis_v2/GEO/validation_data/E_MTAB_386.RData")
colnames(pheno)[1] <- "id"
geneMatrix <- cbind(Reporter.Name=row.names(geneMatrix), geneMatrix)
geneMatrix <- merge(ano[,c("Reporter.Name", "Reporter.Database.Entry.hugo.")], geneMatrix, by = "Reporter.Name")
geneMatrix[1:4, 1:4]
geneMatrix <- geneMatrix[!is.na(geneMatrix$Reporter.Database.Entry.hugo.),]
geneMatrix <- geneMatrix[geneMatrix$Reporter.Database.Entry.hugo. %in% coxGene,]
library(limma)
a <- avereps(as.matrix(geneMatrix[,3:ncol(geneMatrix)]), ID=geneMatrix$Reporter.Database.Entry.hugo.) %>% as.data.frame()
geneMatrix <- cbind(id = geneMatrix[which(rownames(geneMatrix) %in% row.names(a)), 2], a) 

geneMatrix <- t(geneMatrix) %>% as.data.frame()
geneMatrix[1:4, 1:2]
colnames(geneMatrix) <- geneMatrix[1,]
geneMatrix <- geneMatrix[-1,]
geneMatrix[1:4, 1:2]
rownames(geneMatrix) <- gsub("DFCI_", "DFCI-", rownames(geneMatrix))
geneMatrix<- geneMatrix[,] %>% cbind(id=rownames(.), .) %>% as.matrix() %>% merge(pheno, ., by = "id") 
rownames(geneMatrix) <- geneMatrix[,1]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
geneMatrix[1:4, 1:2]
write.table(geneMatrix, "E_MTAB_386_sur_data.txt", col.names = T, sep = "\t")

val <- geneMatrix[,c(2:3, 9:ncol(geneMatrix))] 
a<- row.names(val)

colnames(val)[1:2]<- c("fustime","fustat")

str(val)
val <- apply(val, 2, as.numeric) %>% as.data.frame()
row.names(val) <- a

rs_val=predict(multiCox,type="risk",newdata=val)
risk = as.vector(ifelse(rs_val>median(rs_val),"High_risk","Low_risk"))
val_2 <- cbind(val[,outCol],rs_val,risk) %>% cbind(id=rownames(.), .)
table(val_2$risk)
colnames(val_2)[ncol(val_2)-1] <- "riskScore" 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
write.table(val_2,file="risk.E_MTAB_386.txt", sep="\t",quote=F,row.names=F)

diff=survdiff(Surv(fustime,fustat) ~risk,data = val_2)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(fustime,fustat) ~ risk, data = val_2)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
d <-ggsurvplot(fit,
               pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, 
               risk.table.col = "strata", 
               linetype = "strata", 
               surv.median.line = "hv", 
               ggtheme = theme_bw(),
               palette = c( "#E4392E","#2546F1"),
               break.x.by= 24)
pdf("E_MTAB_386.km.pdf",width=5.5,height=5, onefile = F )
print(d)
dev.off()

#GSE9891
setwd("D:/工作妇肿/课题/disulfidptosis_v2/GEO/validation_data/")
geneMatrix <- read.table("GSE9891 _geneMatrix_ano.txt", sep = "\t", header = T)
pheno <- readxl::read_xlsx("GSE9891.survival.xlsx",1)
XID <- readxl::read_xlsx("GSE9891.ID.xlsx")
colnames(XID)[1:2] <- c("id", "XID")

colnames(pheno)[1] <- "XID"
pheno <- merge(XID[,c(1:2)], pheno, by = "XID")

geneMatrix <- geneMatrix[geneMatrix$gene %in% coxGene,]

geneMatrix <- t(geneMatrix) %>% as.data.frame()
geneMatrix[1:4, 1:2]
colnames(geneMatrix) <- geneMatrix[1,]
geneMatrix <- geneMatrix[-1,]
geneMatrix[1:4, 1:2]
geneMatrix<- geneMatrix %>% cbind(id=rownames(.), .) %>% as.matrix() %>% merge(pheno, ., by = "id") 
rownames(geneMatrix) <- geneMatrix[,1]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
geneMatrix[1:4, 1:3]
write.table(geneMatrix, "GSE9891.txt", col.names = T, sep = "\t")

val <- geneMatrix[,c(9,15, 18:ncol(geneMatrix))] 
colnames(val)[1:2]<- c("fustat", "fustime")

str(val)
id <-rownames(val)
val <- apply(val, 2, as.numeric) %>% as.data.frame()
row.names(val) <- id
outCol=c("fustime","fustat",coxGene)

rs_val=predict(multiCox,type="risk",newdata=val)
risk = as.vector(ifelse(rs_val>median(rs_val),"High_risk","Low_risk"))
val_2 <- cbind(val[,outCol],rs_val,risk) %>% cbind(id=rownames(.), .)
table(val_2$risk)
colnames(val_2)[ncol(val_2)-1] <- "riskScore" 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
write.table(val_2,file="risk.GSE9891.txt", sep="\t",quote=F,row.names=F)

diff=survdiff(Surv(fustime,fustat) ~risk,data = val_2)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(fustime,fustat) ~ risk, data = val_2)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
d <-ggsurvplot(fit,
               pval = TRUE, conf.int = TRUE,
               risk.table = TRUE, 
               risk.table.col = "strata", 
               linetype = "strata",
               surv.median.line = "hv", 
               ggtheme = theme_bw(), 
               palette = c( "#E4392E","#2546F1"),
               xlim= c(0, 100),
               break.x.by= 24)
pdf("GSE9891.km.pdf",width=5.5,height=5, onefile = F )
print(d)
dev.off()


###纳入临床因素做nomogram
rm(list = ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
pheno <- read.table("pheno_colgroup.txt",sep="\t",header=T,check.names=F)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
multiInput=read.table("risk_train.txt", header=T, sep="\t", check.names=F)

trainingdata <- merge(pheno, multiInput, by = "id")
trainingdata[which(trainingdata$stage == 1|trainingdata$stage == 2),4] <- "early"
trainingdata[which(trainingdata$stage == 3|trainingdata$stage == 4),4] <- "late"
trainingdata$HRD <-factor(trainingdata$HRD, levels= c("No", "Yes", "NA"))
trainingdata$BRCAmut <-factor(trainingdata$BRCAmut, levels= c("wt", "Mut", "NA"))
trainingdata$TUMORRESIDUALDISEASE <-factor(trainingdata$TUMORRESIDUALDISEASE, levels= c("R0", "R1", "Not optimal", "NA"))


table(trainingdata$stage)
str(trainingdata)

library(rms)
library (mice)
library(survival)
library(survminer)
library(regplot)

str(pheno)
catVars<-c('stage','TUMORGRADE','TUMORRESIDUALDISEASE', 'ProgressionFreeStatus', 'PlatinumStatus', 'fustat.y')
nonvar<-c('riskScore', 'fustime')
stratavar <- c('risk')
rownames(trainingdata) <-trainingdata$id

trainingdata <- trainingdata[,c('risk', 'riskScore','stage','TUMORGRADE','TUMORRESIDUALDISEASE',
                                'ProgressionFreeStatus', 'PlatinumStatus', 'fustat.y','fustime')]
trainingdata[is.na(trainingdata)] <- "NA"

#Table one
library(tableone)
CreateTableOne(data=trainingdata)
dput(names(trainingdata))


tab2<-CreateTableOne(data=trainingdata, factorVars=catVars, strata = stratavar)
d <- print(tab2,nonnormal = nonvar,showAllLevels = TRUE, noSpaces = TRUE)

write.csv(d, "Table 1.csv")

trainingdata <- trainingdata[,c('riskScore','stage','TUMORGRADE','TUMORRESIDUALDISEASE', 'fustat.y','fustime')]
trainingdata[trainingdata=="NA"] <- NA


##构建基础的cox回归模型

#单因素cox
#install.packages("ezcox")
trainingdata <- na.omit(trainingdata)
str(trainingdata)
trainingdata$TUMORRESIDUALDISEASE <- factor(trainingdata$TUMORRESIDUALDISEASE, levels = c("R0", "R1", "Not optimal"))

library(ezcox)
results <- ezcox(trainingdata, time = "fustime",status = "fustat.y", 
                 covariates = c('riskScore','stage','TUMORGRADE','TUMORRESIDUALDISEASE'))
write.csv(results,file="uniCox_cli.csv",row.names=F,quote=F)

results <- ezcox(trainingdata, time = "fustime",status = "fustat.y", 
                 covariates = c('riskScore','stage','TUMORGRADE','TUMORRESIDUALDISEASE'),
                 return_models=TRUE)
mds = get_models(results)
str(mds, max.level = 1)

pdf("unicox.forest.pdf", width = 8, height = 6)
show_models(mds)
dev.off()


##多因素cox
multiCox_cli=coxph(Surv(fustime, fustat.y) ~riskScore+stage+TUMORGRADE+TUMORRESIDUALDISEASE, data = trainingdata)
multiCox_cli <- step(multiCox_cli,direction = "forward")
multiCoxSum=summary(multiCox_cli)

outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.csv(outTab,file="multiCox_cli.csv",row.names=F,quote=F)

results <- ezcox(trainingdata, time = "fustime",status = "fustat.y", 
                 covariates = c('riskScore'),
                 controls = c('stage','TUMORGRADE','TUMORRESIDUALDISEASE'),
                 return_models=TRUE)
mds = get_models(results)
str(mds, max.level = 1)

pdf("multicox.forest.pdf", width = 8, height = 6)
show_models(mds)
dev.off()

riskScore_cli=predict(multiCox_cli,type="risk",newdata=trainingdata)
trainingdata <- cbind(riskScore_cli=riskScore_cli, trainingdata)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
write.table(trainingdata, "risk.cli_train.txt", sep = "\t")

##new nomogram
table(trainingdata$fustat.y)
library(regplot)
coxfit <- cph(Surv(fustime, fustat.y) ~ riskScore + stage + TUMORGRADE + TUMORRESIDUALDISEASE,
              data = trainingdata, x=T,y=T,surv = T
)

pdf("multicox.cli.nomogram.pdf", width = 8, height = 5)
regplot(coxfit,
        plots = c("density", "boxes"), 
        observation = F, 
        failtime=c(12,36,60), 
        center = T, 
        subticks = T,
        droplines = T,
        title = "nomogram",
        points = T, 
        odds = F, 
        showP = T, 
        rank = "sd", 
        interval="confidence", 
        clickable = F,
        dencol="white",
)
dev.off()

###画训练集的ROC曲线
library(timeROC)
library(survival)
inputFile="risk.cli_train.txt"      
outFile="ROC.cli_train.pdf"         
var="riskScore_cli" 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
multiInput=read.table(inputFile, header=T, sep="\t", check.names=F)
multiInput$fustat <-multiInput$fustat-1
ROC_multiInput=timeROC(T=multiInput$fustime, delta=multiInput$fustat,
                       marker=multiInput[,var], cause=1,
                       weighting='cox',
                       times=c(12,36,60), ROC=TRUE)
pdf(file=outFile,width=5,height=5)
plot(ROC_multiInput,time=12,col='#389613',title=FALSE,lwd=2)
plot(ROC_multiInput,time=36,col='#2546F1',add=TRUE,title=FALSE,lwd=2)
plot(ROC_multiInput,time=60,col='#E4392E',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 12 months: ',sprintf("%.03f",ROC_multiInput$AUC[1])),
         paste0('AUC at 36 months : ',sprintf("%.03f",ROC_multiInput$AUC[2])),
         paste0('AUC at 60 months: ',sprintf("%.03f",ROC_multiInput$AUC[3]))),
       col=c("#389613",'#2546F1','#E4392E'),lwd=2,bty = 'n')
dev.off()

###画验证集的ROC曲线
inputFile="risk.GSE9891.txt"      
outFile="ROC.cli_GSE9891.pdf"         
var="riskScore_cli" 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
geneMatrix <- read.table("GSE9891.txt", header=T, sep="\t", check.names=F)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
multiInput=read.table(inputFile, header=T, sep="\t", check.names=F)

multiInput = merge(geneMatrix,multiInput[,c("id","fustat", "fustime", "riskScore", "risk")], by = "id")
str(multiInput)
multiInput <- multiInput[,c("id", "Stage", "Grade", "ResidualDisease","fustat", "fustime", "riskScore", "risk")]
colnames(multiInput) <- c("id", "stage", "TUMORGRADE", "TUMORRESIDUALDISEASE", "fustat.y", "fustime", "riskScore", "risk")

multiInput[which(multiInput$stage== "I"|multiInput$stage== "II"), 2] <- "early"
multiInput[which(multiInput$stage== "III"|multiInput$stage== "IV"), 2] <- "late"

table(multiInput$TUMORGRADE)
table(trainingdata$TUMORGRADE)
multiInput[which(multiInput$TUMORGRADE== "1"|multiInput$TUMORGRADE== "2"), 3] <- "G2"
multiInput[which(multiInput$TUMORGRADE== "3"), 3] <- "G3"

table(multiInput$TUMORRESIDUALDISEASE)

riskScore_cli=predict(multiCox_cli,type="risk",newdata=multiInput)

multiInput <- cbind(riskScore_cli=riskScore_cli, multiInput)

#绘制
ROC_multiInput=timeROC(T=multiInput$fustime, delta=multiInput$fustat,
                       marker=multiInput[,var], cause=1,
                       weighting='cox',
                       times=c(12,36,60), ROC=TRUE)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
pdf(file=outFile,width=5,height=5)
plot(ROC_multiInput,time=12,col='#389613',title=FALSE,lwd=2)
plot(ROC_multiInput,time=36,col='#2546F1',add=TRUE,title=FALSE,lwd=2)
plot(ROC_multiInput,time=60,col='#E4392E',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 12 months: ',sprintf("%.03f",ROC_multiInput$AUC[1])),
         paste0('AUC at 36 months : ',sprintf("%.03f",ROC_multiInput$AUC[2])),
         paste0('AUC at 60 months: ',sprintf("%.03f",ROC_multiInput$AUC[3]))),
       col=c("#389613",'#2546F1','#E4392E'),lwd=2,bty = 'n')
dev.off()


###画验证集E_MTAB_386的ROC曲线
inputFile="risk.E_MTAB_386.txt"      
outFile="ROC.cli_E_MTAB_386.pdf"         
var="riskScore_cli" 

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GEO/validation_data/")
load("D:/工作妇肿/课题/disulfidptosis_v2/GEO/validation_data/E_MTAB_386.RData")
geneMatrix <- pheno
colnames(geneMatrix)[1] <- "id"

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/validation/")
multiInput=read.table(inputFile, header=T, sep="\t", check.names=F)

multiInput = merge(geneMatrix,multiInput[,c("id", "riskScore", "risk")], by = "id")
colnames(multiInput)
multiInput <- multiInput[,c(1:3, 6:10)]
colnames(multiInput) <- c("id", "fustime", "fustat.y", "stage", "TUMORRESIDUALDISEASE",  "TUMORGRADE", "riskScore", "risk")
str(multiInput)

multiInput$fustime <- as.numeric(multiInput$fustime)

table(multiInput$stage)
multiInput$stage <- as.character(multiInput$stage)
multiInput[which(multiInput$stage=="Iic"),4] <- "early"
multiInput[which(multiInput$stage=="IIIa"),4] <- "late"
multiInput[which(multiInput$stage=="IIIb"),4] <- "late"
multiInput[which(multiInput$stage=="IIIc"),4] <- "late"
multiInput[which(multiInput$stage=="IV"),4] <- "late"

table(multiInput$TUMORRESIDUALDISEASE)
multiInput$TUMORRESIDUALDISEASE <- as.character(multiInput$TUMORRESIDUALDISEASE)
multiInput[which(multiInput$TUMORRESIDUALDISEASE=="No"),5] <- "Not optimal"
multiInput[which(multiInput$TUMORRESIDUALDISEASE=="Yes"),5] <- "R0"
multiInput[which(multiInput$TUMORRESIDUALDISEASE=="Unknown"),5] <- NA

table(multiInput$TUMORGRADE)
multiInput[which(multiInput$TUMORGRADE== "high"), 6] <- "G3"

multiInput <- na.omit(multiInput)

riskScore_cli=predict(multiCox_cli,type="risk",newdata=multiInput)

multiInput <- cbind(riskScore_cli=riskScore_cli, multiInput)

#绘制
ROC_multiInput=timeROC(T=multiInput$fustime, delta=multiInput$fustat,
                       marker=multiInput[,var], cause=1,
                       weighting='cox',
                       times=c(12,36,60), ROC=TRUE)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
pdf(file=outFile,width=5,height=5)
plot(ROC_multiInput,time=12,col='#389613',title=FALSE,lwd=2)
plot(ROC_multiInput,time=36,col='#2546F1',add=TRUE,title=FALSE,lwd=2)
plot(ROC_multiInput,time=60,col='#E4392E',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 12 months: ',sprintf("%.03f",ROC_multiInput$AUC[1])),
         paste0('AUC at 36 months : ',sprintf("%.03f",ROC_multiInput$AUC[2])),
         paste0('AUC at 60 months: ',sprintf("%.03f",ROC_multiInput$AUC[3]))),
       col=c("#389613",'#2546F1','#E4392E'),lwd=2,bty = 'n')
dev.off()


###算死亡评分
rm(list = ls())  
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
exp <- read.table(file ="log2cpm___new.txt", header=T,check.names=F, sep = "\t")

head(exp[1:4,1:4])
rt <- exp
rt[1:4, 1:4]

row.names(rt) <- rt[,1]
geneList <- rt[,1]

geneList <- sapply(strsplit(geneList, "\\."), "[", 1)

library(clusterProfiler)

gene.df <- bitr(geneList, fromType="ENSEMBL",
                toType="SYMBOL", 
                OrgDb = "org.Hs.eg.db")
gene.df <- na.omit(gene.df)

rt2 <- rt %>% as.data.frame() %>% cbind(ENSEMBL=rownames(.), .)

rt2[1:4, 1:4]
rt2$ENSEMBL <- sapply(strsplit(rt2$ENSEMBL, "\\."), "[", 1)

rt3<- merge(gene.df, rt2, by= "ENSEMBL")
rt3[1:4, 1:4]
rt3<- as.matrix(rt3)
row.names(rt3)<- rt3[,2]
rt3<- rt3[,-c(1:3)]
rt3[1:4, 1:4]

library(limma)
library(magrittr)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")

#对同一基因取均值去重
rt4<-avereps(as.matrix(rt3))
head(rt4[1:4,1:4])
dim(rt4)
total_gene <- rownames(rt4)

mer3 <- rt4 %>% t() %>% cbind(id=row.names(.), .)
mer3[1:4, 1:4]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
multiInput=read.table("risk_train.txt", header=T, sep="\t", check.names=F)
multiInput<- multiInput[,-c(4:(ncol(multiInput)-2))]
mer4 <- merge(multiInput, mer3, by="id")
mer4[1:4,1:4]
dim(mer4)

#画4个关键基因的组间对比
library(reshape2)
mer4 <- as.data.frame(mer4)
four_genes <- mer4[,c("risk", "riskScore", "SLC7A11", "SLC3A2", "RPN1", "NCKAP1")]
four_genes_2 <- melt(four_genes, id.vars = c("risk", "riskScore"), variable.name = "gene", value.name = "expr")
four_genes_2$expr <- as.numeric(four_genes_2$expr)

pdf(file="Disulfidptosis_fourgenes_boxplot.pdf", width = 5.5, height = 4.5)
ggplot(four_genes_2, mapping = aes(x = gene, y = expr, color=factor(risk)))+
  geom_boxplot(aes(x = gene, y = expr, color=factor(risk)),
               stat = "boxplot",
               alpha = 1, 
               width = 1,
               outlier.color = "grey",
               outlier.fill = "grey", 
               outlier.shape = 19, 
               outlier.size = 1, 
               outlier.alpha = 0.1) + 
  stat_compare_means(label = "p.signif", label.x=1.3, hide.ns = T)+ 
  scale_fill_manual(values = c("Low_risk"= "#2546F1", "High_risk"="#E4392E"))+
  scale_color_manual(values = c("Low_risk"= "#2546F1", "High_risk"="#E4392E"))+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=30, hjust = 1, vjust=1))
dev.off()


row.names(mer4) <- mer4$id
exp<- mer4[,-c(1:5)]
exp[1:5, 1:5]
rs<- mer4[,c(1:5)]
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
write.table(rs,file="coxriskscore.txt",sep="\t",quote=F,row.names=T)
write.table(exp, file="TCGA_ov_exp_log2cpm.txt",sep="\t",quote=F,row.names=T)

##计算双硫死亡的得分表
exp[1:4,1:4]
b=as.matrix(t(exp))
dim(b)
a <- apply(b, 2, as.numeric) %>% as.data.frame() #转化为numeric格式
row.names(a) <- row.names(b)
colnames(a) <- colnames(b)
dim(a)
a[1:4, 1:4]

library(GSVA)
library(dplyr)
genes<-read.table("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/disulfidptosis_score.txt")
setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/")

gs<-data.frame(na.omit(genes$V1))
gs=as.list(gs)
gs=lapply(gs, function(x)x[!is.na(x)])
ssgsea_score=gsva(as.matrix(a),gs, method="ssgsea",ssgsea.norm=T,verbose=T)
file_name <- "disulfidptosis_ssGSEA.csv"
write.csv(ssgsea_score,file_name)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/")
rt<-read.csv("disulfidptosis_ssGSEA.csv", check.names = F)
rt<- t(rt)
colnames(rt)<- rt[1,]
rt<- rt[-1,]
rt2<- as.data.frame(rt)
colnames(rt2)[1] <- "disulfidptosis_score"
rt2 <- sapply(rt2,as.numeric,2)
rt2[1:4,]

ds<-ifelse(rt2[,"disulfidptosis_score"]<median(rt2[,"disulfidptosis_score"]),"Low_score","High_score")

t1<-cbind(rt,ds_class=ds)
t1<-cbind(id=row.names(t1),t1)
colnames(t1)[2] <- "ds"

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
rs<- read.table("coxriskscore.txt",header=T,sep="\t",check.names=F)
t2<- merge(rs,t1,by="id")
setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/")
write.table(t2,file = "Tumor_DisulfidptosisScoreGroup_Survival.txt",sep="\t",row.names = F)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/")
t2<- read.table("Tumor_DisulfidptosisScoreGroup_Survival.txt",header=T,sep="\t",check.names=F) ##包含cox风险分层、5种死亡评分及分层

library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(ggpubr)

t2$risk <- as.factor(t2$risk)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSVA/")
pdf(file="DisulfidptosisScore_boxplot.pdf", width = 5.5, height = 4.5)
ggplot(t2, mapping = aes(x = factor(risk), y = ds, colour=factor(risk))) + 
  geom_quasirandom( aes(x = factor(risk), y = ds, fill = risk,colour=factor(risk)),
                    alpha = 1, 
                    width = 0.4, 
                    size = 0.8)+ 
  stat_boxplot(geom = "errorbar", 
               width = 0.1, 
               size = 0.1, 
               color = "black")+
  geom_boxplot(aes(x = factor(risk), y = ds),
               stat = "boxplot",
               alpha = 1, 
               width = 0.1, 
               color = "black",
               outlier.color = "grey", 
               outlier.fill = "grey", 
               outlier.shape = 19, 
               outlier.size = 1, 
               outlier.alpha = 0.1) + 
  stat_compare_means(label = "p.signif", label.x=1.3, hide.ns = T)+ 
  scale_fill_manual(values = c("Low_risk"= "#2546F1", "High_risk"="#E4392E"))+
  scale_color_manual(values = c("Low_risk"= "#2546F1", "High_risk"="#E4392E"))+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=30, hjust = 1, vjust=1))
dev.off()

#差异基因分析
rm(list = ls())
library(magrittr)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
exp <- read.table("TCGA_ov_exp_log2count.txt",sep="\t",check.names = F,header = T)
exp[1:5,1:5]
str(exp)

exp1 <- rbind(colnames(exp), exp)
head(exp1[1:4,1:4])

exp2<- t(exp1)
head(exp2[1:4,1:4])
exp2 <- exp2[,] %>% as.data.frame(.,  stringsAsFactors = F)
head(exp2[1:4,1:4])
colnames(exp2)[1] <- "id"
str(exp2)


library(limma)
library(magrittr)
exp3 <-avereps(as.matrix(exp2[,-1]), ID=exp2$id)
exp3 <- as.data.frame(exp3, stringsAsFactors = F)
dim(exp3)
str(exp3)

exp3[1:4, 1:4] #exp3还是log2(count+1)格式
gene <- row.names(exp3)
sampleID <- colnames(exp3)
exp4 <- apply(exp3, 2, as.numeric) %>% as.data.frame()
exp4[1:4, 1:4]

row.names(exp4) <- gene

exp5 <- 2^exp4-1 #exp5 转为count格式

exp5 <- exp5 %>% t() %>% data.frame()
exp5[1:4, 1:4]

row.names(exp5) <- sampleID
exp5 <- cbind(id= rownames(exp5),exp5)
exp5[1:4, 1:4]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
anno2<- read.table("coxriskscore.txt",header=T,sep="\t",check.names=F)

exp_cluster <- merge(anno2,exp5, by = "id")
exp_cluster[1:4, 1:4]
exp_cluster <- exp_cluster[order(exp_cluster$risk, decreasing = F), ]
row.names(exp_cluster) <- exp_cluster$id

exp <- t(exp_cluster[,6:ncol(exp_cluster)])
gene <- row.names(exp)
exp <- apply(exp, 2, as.numeric) %>% as.data.frame()
row.names(exp) <- gene
exp[1:4, 1:4]

exp <- round(exp)

exp <- exp[apply(exp, 1, function(x) sum(x > 1) > 45*0.75), ]
dim(exp)
write.table(exp, "exp_cluster_count_new.txt", sep="\t",row.names = T, quote=F)

exp_anno <- exp_cluster[,c(1,5)]

library(limma)
library(edgeR)

#分组矩阵design构建
design <- model.matrix(~0+factor(exp_anno$risk)) 
table(exp_anno$risk)

colnames(design) <- c("High_risk", "Low_risk")
rownames(design) <- colnames(exp)

## 表达矩阵DGEList构建与过滤低表达基因
dge <- DGEList(counts=exp) 
keep.exprs <- filterByExpr(dge,design=design) 
table(keep.exprs)
dge <- dge[keep.exprs,keep.lib.sizes=FALSE] 
dge <- calcNormFactors(dge)  
cont.matrix <- makeContrasts(High_risk-Low_risk, 
                             levels = design)

de <- voom(dge,design,plot=TRUE, normalize="quantile")
fit1 <- lmFit(de, design)              
fit2 <- contrasts.fit(fit1,cont.matrix) 
efit <- eBayes(fit2, trend=F)  #Apply empirical Bayes smoothing to the standard errors

tempDEG <- topTable(efit, n=Inf, adjust= "fdr") 
DEG_limma_voom  <- na.omit(tempDEG)
head(DEG_limma_voom)

diffsig <- DEG_limma_voom 
setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
write.csv(diffsig, "all.limmaOut.csv")

foldChange = 0.5
padj = 0.05
All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
dim(All_diffSig)
write.csv(All_diffSig, "all.diffsig_p0.05.csv") 

diffup <-  All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup.csv")
#
diffdown <- All_diffSig[(All_diffSig$adj.P.Val < padj & (All_diffSig$logFC < -foldChange)),]
write.csv(diffdown, "diffdown.csv")


rm(list = ls())
## 提取差异基因的表达量
setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
load("TCGA_OV_count_transformation.Rdata")
exp<- vsdd

DEG_id <- read.csv("all.diffsig_p0.05.csv", header = T) 
head(DEG_id)
DEG_id2 <- DEG_id[order(abs(DEG_id$logFC), decreasing = T),]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
pheno6 <- read.table("BRCA_hrd_tcga.txt",sep = "\t", check.names = F)
pheno6[is.na(pheno6)] <- "NA"

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
anno2<- read.table("coxriskscore.txt",header=T,sep="\t",check.names=F)

DEG_id <- unique(DEG_id2$X)
hmexp <- na.omit(exp)
hmexp[1:4, 1:4]
hmexp <- hmexp[DEG_id[1:40],]
str(hmexp)

library(pheatmap)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
colgroup <- read.table("pheno_colgroup.txt",sep="\t",header=T,check.names=F)


setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
## 样本注释信息 
annotation_col <- data.frame(colgroup[,c(2:5, 7:ncol(colgroup))])
str(annotation_col)

table(annotation_col$fustat)
annotation_col$fustat <- as.factor(annotation_col$fustat)
rownames(annotation_col) <- colnames(hmexp)

col = colorRampPalette(c("#1435AD",  "white", "#FFB100"))(300)

table(annotation_col$cluster)
table(annotation_col$ProgressionFreeStatus)

annotation_col$age <- as.numeric(annotation_col$age)
ann_colors = list( cluster = c(High_risk = "#F34800", Low_risk= "#00447E"), 
                   stage = c(`1`="#e9f679", `2`= "#9bdf46", `3`= "#25a55f", `4`= "#346473", `NA` = "#eaeaea"),
                   age = c("white", "firebrick"),
                   HRD = c(No = "#08d9d6", Yes = "#ff2e63", `NA` = "#eaeaea"),
                   BRCAmut = c(Mut = "#ff6f3c", wt = "#155263", `NA` = "#eaeaea"),
                   fustat = c(`2` = "#476269", `1` = "#f5e1da"),
                   TUMORGRADE = c(G2 = "#a7d7c5", G3 = "#5c8d89", `NA` = "#eaeaea"),
                   TUMORRESIDUALDISEASE = c(R0="#fbd0f5", R1="#f7f680", `Not optimal`= "#94f6f2", `NA` = "#eaeaea"),
                   ProgressionFreeStatus = c(DiseaseFree="#576aff", `Recurred/Progressed` = "#ffc688"),
                   PlatinumStatus = c(Resistant="#edff91", Sensitive="#5ac7a2", Missing="#eaeaea"))

##  绘制热图 
setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
bk = unique(c(seq(-3,3, length=300)))
pdf(file = "deg_heatmap.pdf", height = 10, width = 15)
pheatmap(hmexp,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         breaks = bk,
         color = col,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         scale = "row", ## none, row, column
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 6,
         gaps_col = F,
         border = FALSE)
dev.off()

###画valcano图
## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
diffsig <- read.csv("all.limmaOut.csv", header = T)
row.names(diffsig) <- diffsig[,1]
diffsig <- diffsig[,-1]

# 开始绘图
foldChange = 0.5
padj = 0.05
label = subset(diffsig, adj.P.Val < padj & abs(logFC)>foldChange)
label1 = rownames(label)

Significant=ifelse((diffsig$adj.P.Val < padj & abs(diffsig$logFC) > foldChange), ifelse(diffsig$logFC > foldChange,"Up","Down"), "Not")

P <- ggplot(diffsig, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept = c(-foldChange, foldChange),colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(padj),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig)+
  xlim(-2.5,2.5)+ylim(0,25)+theme_bw()

diffsig$label=ifelse(diffsig$adj.P.Val < padj & abs(diffsig$logFC)>foldChange,rownames(diffsig),"")

pdf('volcano.pdf',width = 7,height = 6.5)  
P+geom_text_repel(data = diffsig, aes(x = logFC, 
                                      y = -log10(adj.P.Val), 
                                      label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
dev.off()

###分析临床因素
rm(list = ls())
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(ggpubr)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/TCGA_ov/")
pheno <- read.table("pheno_colgroup.txt",sep="\t",header=T,check.names=F)
pheno <- pheno[,-2]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
risk <- read.table("coxriskscore.txt",sep="\t",header=T,check.names=F)
pheno <- merge(pheno, risk[,c(1,4,5)], by = "id")

pheno$age <- as.numeric(pheno$age)
pheno$TUMORRESIDUALDISEASE <- factor(pheno$TUMORRESIDUALDISEASE, levels = c("R0", "R1", "Not optimal"))
pheno$fustat <- factor(pheno$fustat, levels = c("2", "1"))
pheno$BRCAmut <- factor(pheno$BRCAmut, levels = c("Mut", "wt"))
pheno$HRD <- factor(pheno$HRD, levels = c("Yes", "No"))

colnames(pheno)

Age <- ggplot(data=pheno, mapping=aes(x=as.factor(risk),y=age, fill= as.factor(risk)))+
  geom_violin(color= "black", width=0.5, alpha = 1,trim = F)+ 
  geom_boxplot(stat = "boxplot",
               alpha = 0.5, 
               width = 0.1, 
               fill = "white",
               outlier.color = "grey", 
               outlier.fill = "grey", 
               outlier.shape = 19, 
               outlier.size = 1, 
               outlier.alpha = 0.1) + 
  stat_compare_means(label = "p.signif", label.x=1.3, hide.ns = F)+ 
  scale_fill_manual(values = c("Low_risk"= "#709AE1FF", "High_risk"="#FED439FF"))+   
  scale_color_manual(values = c("Low_risk"= "#709AE1FF", "High_risk"="#FED439FF"))+
  theme_bw()+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=30, hjust = 1, vjust=1))
Age
dev.off()

setwd("D:/工作妇肿/课题/disulfidptosis_v2/clinical_relavance/")
pdf("age_violin.pdf", width = 3, height = 4.5)
Age
dev.off()

library(grid)
library(dplyr)
colnames(pheno)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12)
)


setwd("D:/工作妇肿/课题/disulfidptosis_v2/clinical_relavance/")
vars1 <- c("stage", "TUMORRESIDUALDISEASE","TUMORGRADE", "fustat", "HRD", "BRCAmut")
table(pheno$TUMORRESIDUALDISEASE)
for (i in vars1){
  plot <- pheno[which(!is.na(pheno[,i])),]
  n <- plot[,c("risk",i)]
  n[,i] <- as.factor(n[,i])
  colnames(n)[2] <- "var"
  m <- paste0(i, ".pdf")
  i_surbar <- ggplot(n, aes(x = risk, group = var, fill = var))+
    geom_bar(position="fill", stat = "count") +
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = c("#FED439FF","#709AE1FF", "#FD7446FF", "#D2AF81FF"))+theme_bw()+mytheme
  
  y <- paste0(i, "_surbar.pdf")
  i_surbar;ggsave(y,i_surbar,height=8,width=5,unit="cm")
}

extrafont::loadfonts()

library(ggstatsplot)
library(hrbrthemes)
for (i in vars1){
  plot <- pheno[,c("risk",i)]
  n <- plot[which(!is.na(pheno[,i])),]
  n[,i] <- as.factor(n[,i])
  colnames(n)[2] <- "var"
  m <- paste0(i, ".pdf")
  i_surbar <- ggbarstats(
    data = n,
    x = var,
    y = risk,
    title = i,
    xlab = "Cluster",
    legend.title = i,
    #ggtheme = hrbrthemes::theme_ipsum(base_family = "Roboto Condensed"),
    ggplot.component = list(ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))),
    bf.message = FALSE
  )+ scale_fill_manual(values = c("#FED439FF","#709AE1FF", "#FD7446FF", "#D2AF81FF"))
  y <- paste0(i, "_statbar.pdf")
  pdf(y,height=4,width=2.5)
  print(i_surbar)
  dev.off()
}

#TCGA_mut
library(TCGAbiolinks)
a <- GDCquery(project = "TCGA-OV",
              data.category = "Copy Number Variation",
              data.type = "Masked Copy Number Segment",
              workflow.type = "DNAcopy")
setwd("G:/ov/")
GDCdownload(a, method = "api")

cnv_files <-GDCprepare(query = a, save = TRUE, save.filename = "CNV_TCGA_OV.rda")

cnv<- load("G:/ov/CNV_TCGA_OV.rda")
hnsc_seg <- eval(parse(text=cnv))
hnsc_seg <- hnsc_seg[,-1]
hnsc_seg <- hnsc_seg[,c('Sample', 'Chromosome', 'Start', 'End', "Num_Probes", "Segment_Mean")]

tumor_seg <- hnsc_seg[substr(hnsc_seg$Sample,14,15) == "01",]
write.table(tumor_seg, "OV_tumor_seg.txt", sep = "\t", col.names = TRUE, row.names = F)

hg_marker_file <- read.delim("snp6.na35.remap.hg38.subset.txt.gz")
hg_marker_file <- hg_marker_file[hg_marker_file$freqcnv == "FALSE", ]
hg_marker_file <- hg_marker_file[,c(1,2,3)]
write.table(hg_marker_file, "hg_marker_file.txt", sep = "\t", col.names = TRUE, row.names = F)

setwd("G:/ov/")
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

snp.data <- GDCprepare(query, save = T,save.filename = "TCGA-OV_SNP.rda")


library(maftools)

load(file = "G:/ov/TCGA-OV_SNP.rda")

maf.ov <- data 
maf.ov
maf <- read.maf(maf.ov)

getSampleSummary(maf)

setwd("G:/ov/")
#Read GISTIC results along with MAF
colnames(risk)[1] <- "Tumor_Sample_Barcode"
ov.plus.gistic = read.maf(
  maf = maf.ov,
  gisticAllLesionsFile = 'G:/ov/all_leisions.conf_90.txt',
  gisticAmpGenesFile = 'G:/ov/amp_genes.conf_90.txt',
  gisticDelGenesFile = 'G:/ov/del_genes.conf_90.txt',
  gisticScoresFile = 'G:/ov/scores.gistic',
  isTCGA = TRUE,
  verbose = FALSE, 
  clinicalData = risk
)
oncoplot(maf = ov.plus.gistic, top = 10)
save(ov.plus.gistic, file="ov.plus.gistic.RData")

rm(list = ls())
setwd("G:/ov/")
load("ov.plus.gistic.RData")  
risk <- read.table("disulfidptosis/coxriskscore.txt")
genes <- read.table("disulfidptosis/disulfidptosis.txt")
genes = genes$V1

library(maftools)
risk_mut <- subsetMaf(ov.plus.gistic, tsb = risk$id)
high_risk <- subsetMaf(ov.plus.gistic, tsb = risk$id[which(risk$risk=="High_risk")])
low_risk <- subsetMaf(ov.plus.gistic, tsb = risk$id[which(risk$risk=="Low_risk")])
oncoplot(maf = high_risk, top = 10)
oncoplot(maf = low_risk, top = 10)

oncoplot(maf = high_risk, genes = genes)
oncoplot(maf = low_risk, genes = genes)

#添加临床信息注释
getClinicalData(x = risk_mut)

fabcolors = c("red", "blue")
names(fabcolors) = c("High_risk", "Low_risk")
fabcolors = list(risk = fabcolors)

print(fabcolors)

oncoplot(maf = risk_mut, top = 20, clinicalFeatures = 'risk')

pdf("oncoplot_risk_cnv.snp.pdf", width = 12, height = 6)
oncoplot(
  maf = risk_mut, 
  top = 20, 
  clinicalFeatures = 'risk',
  sortByAnnotation = TRUE,
  annotationColor = fabcolors
)
dev.off()


##gsea分析
rm(list = ls())

library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
DEG_id <- read.csv("all.limmaOut.csv", header = T)  
head(DEG_id)

DEG_id <- DEG_id[order(DEG_id$logFC, decreasing = T),]
gene_list <- DEG_id$logFC
names(gene_list) <- DEG_id$X
print(head(gene_list))

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
anno2<- read.table("coxriskscore.txt",header=T,sep="\t",check.names=F)

gene_list<-gene_list[!is.na(names(gene_list))]
gene_list_up <- gene_list[gene_list>0]
gene_list_down <- gene_list[gene_list<0]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
gene_set<-read.gmt("h.all.v2022.1.Hs.symbols.gmt")
print(gene_set)

gsea.raw<-GSEA(geneList = gene_list,
               nPerm= 1000,
               pvalueCutoff = 1,
               TERM2GENE=gene_set)
print(head(gsea.raw@result))


gsea <- gsea.raw@result
gsea2 <- gsea[abs(gsea$NES)>1 & gsea$pvalue < 0.05,]
colnames(gsea2)

down_gsea <- gsea2 %>%
  filter(NES < -1) %>% 
  filter(pvalue<0.05) %>% 
  mutate(group=-1) 

up_gsea <- gsea2 %>%
  filter(NES > 1) %>% 
  filter(pvalue<0.05) %>% 
  mutate(group=1) 
gsea2 <- rbind(up_gsea, down_gsea)

gsea2$pvalue = -log10(gsea2$pvalue)
gsea2$NES=gsea2$pvalue*gsea2$group 
gsea2=gsea2[order(gsea2$NES,decreasing = T),]

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
save(gsea2, gsea.raw, file="gsea_hallmark.RData")

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
load("gsea_hallmark.RData")
gsea2$Description <- gsub("HALLMARK_", "", gsea2$Description)
g<- ggplot(gsea2, aes(x=reorder(Description,order(NES, decreasing = T)), y=NES, fill=group)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="#709AE1FF",high="#FED439FF",guide = FALSE) +
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="NES") +
  coord_flip() +
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Hallmark Pathway Enrichment") 

pdf("hallmark_gsea.pdf", width = 10, height = 6)
g
dev.off()

geneSet <- c("HALLMARK_HYPOXIA",
             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
             "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
             "HALLMARK_ANGIOGENESIS", 
             "HALLMARK_TGF_BETA_SIGNALING")

gsea.raw@result$ID
geneSet_ID <- which(gsea.raw@result$ID %in% geneSet)

library(ggsci)
library(enrichplot)
col_gsea1<-pal_simpsons()(5)
gseaplot<- gseaplot2(gsea.raw, geneSetID = geneSet_ID,
                     title = "",
                     color = col_gsea1[1:5],
                     base_size = 14,
                     rel_heights = c(1, 0.2, 0.2),
                     subplots = 1:3,
                     pvalue_table = T,
                     ES_geom = "line"
)

pdf("hallmark_gseaplot.pdf", width = 10, height = 8)
gseaplot
dev.off()

#low_risk组高表达的通路
geneSet <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
             "HALLMARK_E2F_TARGETS", 
             "HALLMARK_MYC_TARGETS_V1", 
             "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
             "HALLMARK_INTERFERON_GAMMA_RESPONSE")

gsea.raw@result$ID
geneSet_ID <- which(gsea.raw@result$ID %in% geneSet)

library(ggsci)
library(enrichplot)
col_gsea1<-pal_simpsons()(5)
gseaplot<- gseaplot2(gsea.raw, geneSetID = geneSet_ID,
                     title = "",
                     color = col_gsea1[1:5],
                     base_size = 14,
                     rel_heights = c(1, 0.2, 0.2),
                     subplots = 1:3,
                     pvalue_table = T,
                     ES_geom = "line"#line or dot
)

pdf("hallmark_gseaplot_lowrisk.pdf", width = 10, height = 8)
gseaplot
dev.off()


#actin related pathways
actin_pathways <- c("GOMF_ACTIN_BINDING", "GOMF_CADHERIN_BINDING", "GOMF_ACTIN_FILAMENT_BINDING",
                    "GOMF_CELL_ADHESION_MOLECULE_BINDING", 	"GOMF_CYTOSKELETAL_PROTEIN_BINDING" )

gsea3 <- gsea[which(gsea$Description %in% actin_pathways),]

gsea.raw@result <- gsea3

library(ggsci)
col_gsea1<-pal_simpsons()(3)
gseaplot<- gseaplot2(gsea.raw, geneSetID = gsea3$Description,
                     title = "",
                     color = col_gsea1[1:3],
                     base_size = 14,
                     rel_heights = c(1, 0.2, 0.2),
                     subplots = 1:3,
                     pvalue_table = T,
                     ES_geom = "line"#line or dot
)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
pdf("actin_related_pathways_gseaplot.pdf", width = 10, height = 8)
gseaplot
dev.off()


library(clusterProfiler)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/DEG/")
diffsig <- read.csv("all.diffsig_p0.05.csv", header = T)  
row.names(diffsig) <- diffsig[,1]
diffsig <- diffsig[,-1]

foldChange = 0.5
padj = 0.05

gene_up = rownames(diffsig[diffsig$logFC > foldChange & diffsig$adj.P.Val< padj,]) 
gene_down = rownames(diffsig[diffsig$logFC < foldChange & diffsig$adj.P.Val< padj,]) 
gene_list = c(gene_up,gene_down)

gene.df <- bitr(gene_list, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
gene.df <- na.omit(gene.df)

gene_up <- gene.df[which(gene.df$SYMBOL %in% gene_up),]
gene_down <- gene.df[which(gene.df$SYMBOL %in% gene_down),]

GO <- enrichGO(gene = gene.df$ENTREZID, 
               OrgDb = "org.Hs.eg.db", 
               pvalueCutoff = 0.05,
               ont="all")

a <- GO@result
setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
pdf("GO_dotplot.pdf")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

#KEGG
KEGG_database <- 'hsa' 
KEGG<-enrichKEGG(gene.df$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.5,
                 minGSSize = 10, maxGSSize = 5000
)

ekegg <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_diff_dt <- data.frame(ekegg)
p1 <- dotplot(ekegg, showCategory=20, title = 'KEGG Pathway')

pdf("KEGG_dotplot.pdf", width = 8, height = 4)
p1
dev.off()

save(GO, KEGG, file="GO_KEGG.RData")



#免疫特征分析
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
train<-read.table("coxriskscore.txt",header=T,sep="\t")

setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/timer/")
timer <- read.csv("infiltration_estimation_for_tcga.csv")
a <- substr(timer$cell_type, 14, 15)
table(a)
timer <- cbind(a, timer)
timer <- timer[!a==11,]
timer$cell_type <- substr(timer$cell_type, 1, 12)
timer <- timer[which(timer$cell_type %in% train$id),]
table(timer$a)

timer2 <- timer[!timer$a=="02",]
row.names(timer2) <- timer2$cell_type
timer2 <- timer2[,-1]
colnames(timer2)[1] <- "id"

timer3 <- merge(train[,c(1,4,5)], timer2, by = "id")

library(psych)
res = corr.test(timer3[,2], timer3[,4:ncol(timer3)],method = "pearson", adjust = "fdr")
a <- res$r
b <- res$p.adj
c <- res$ci2
d<- rbind(a,b)
row.names(d) <- c("ce", "p.adj")
d <- t(d) %>% as.data.frame()
d <- d[d$p.adj<0.05,] %>% na.omit()
d <- cbind(cell_type=row.names(d),d)
cell_type <- sapply(strsplit(d$cell_type, "\\_"), "[", 1)
algorithm <- sapply(strsplit(d$cell_type, "\\_"), "[", 2)
d <- cbind(new_celltype=cell_type, algorithm= algorithm, d)

library(dplyr)

plot_order <- d %>% 
  group_by(new_celltype) %>% 
  summarise(m = max(ce)) %>% 
  arrange(desc(m)) %>% 
  pull(new_celltype)

## `summarise()` ungrouping output (override with `.groups` argument)
d$new_celltype <-  factor(d$new_celltype, levels = plot_order)
a <- -log10(d$p.adj)
for(i in 1:length(a)){
  if(a[i]>3){
    a[i] <- "***"
  }else if(a[i] > 2 ){
    a[i] <- "**"
  }else if(a[i] > -log10(0.05)){
    a[i] <- "*"
  }else (a[i] <- "")
}
d$sig <- a

library(ggplot2)
library(ggpubr)
pdf("sig_celltype.pdf", width = 5, height = 6)
ggdotchart(d, x = "new_celltype", y = "ce",
           color = "algorithm",                                
           palette = c("#00AFBB", "#E7B800", "#FC4E07", "#2546F1", "#E4392E", "#389613", "#AA6C57"), 
           sorting = "descending",                      
           add = "segments",                             
           add.params = list(color = "lightgray", size = 1.5),
           rotate = TRUE,                              
           #group = "algorithm",                                
           dot.size = 6,                                
           label = d$sig,                     
           font.label = list(color = "white", size = 8, 
                             vjust = 0.5),               
           ggtheme = theme_pubr(),                       
           xlab=""
)+
  geom_hline(aes(yintercept=0), color="lightgray")

timer4 <- timer3[,c("id", "risk", d$cell_type)]
colnames(timer4)[1] <- "sample"
library(reshape2)
TME_New = melt(data=timer4,id=c("sample","risk"),measure.vars=colnames(timer4)[3:ncol(timer4)], variable.name = "Celltype", value.name = "Composition")
head(TME_New)

TME_New$Celltype <- as.character(TME_New$Celltype)
cell_type <- sapply(strsplit(TME_New$Celltype, "\\_"), "[", 1)
algorithm <- sapply(strsplit(TME_New$Celltype, "\\_"), "[", 2)
TME_New <- cbind(new_celltype=cell_type, algorithm= algorithm, TME_New)
unique(cell_type)

lymphoid <- c("T.cell.CD4.", "T.cell.CD8.", "T.cell", "T.cell.CD4..memory", "T.cell.NK", "T.cell.CD4..Th1", "NK.cell")
myeloid <- c("Macrophage", "Myeloid.dendritic.cell.activated", "Neutrophil", "Monocyte", "Macrophage.Monocyte")
stem_cell <- c("Common.lymphoid.progenitor", "Common.myeloid.progenitor", "Hematopoietic.stem.cell")
stromal_cell <- c("Cancer.associated.fibroblast", "Endothelial.cell")
others <- c("stroma.score", "microenvironment.score", "uncharacterized.cell")

for(i in 1:nrow(TME_New)){
  if (TME_New$new_celltype[i] %in% lymphoid){
    TME_New$class[i] <- "lymphoid"
  } else if(TME_New$new_celltype[i] %in% myeloid){
    TME_New$class[i] <- "myeloid"
  } else if(TME_New$new_celltype[i] %in% stem_cell){
    TME_New$class[i] <- "stem_cell"
  } else if(TME_New$new_celltype[i] %in% stromal_cell){
    TME_New$class[i] <- "stromal_cell"
  } else if(TME_New$new_celltype[i] %in% others){
    TME_New$class[i] <- "others"
  } else TME_New$class[i] <- NA
}

TME_New$Celltype <- factor(TME_New$Celltype)
TME_New$Composition <- as.numeric(TME_New$Composition)
TME_New$risk <- factor(TME_New$risk, levels = c("Low_risk", "High_risk"))
TME_New$class <- factor(TME_New$class, levels = c("others", "stem_cell", "stromal_cell", "myeloid", "lymphoid"))
TME_New <- TME_New[order(TME_New$class, decreasing = F),]

library(Rmisc)
TME_summary <- summarySE(TME_New, measurevar="Composition", groupvars=c("risk","Celltype"))
head(TME_summary)

library(ggplot2)
library(ggpubr)
library(ggunchained)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/timer/")
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }


library(forcats) ##使用fct_reorder()给图进行排序
box_TME1 <- TME_New[!TME_New$algorithm == "MCPCOUNTER",] %>%
  mutate(Celltype = factor(Celltype, levels = unique(TME_New$Celltype)))%>% 
  ggplot(aes(x = Celltype, y = Composition, color = risk))+ 
  labs(y="Tumor Environment Score",x= NULL,title = NULL)+  
  geom_boxplot(position=position_dodge(0.5),width=0.5,alpha = 0.8, outlier.size = 0)+ 
  scale_color_manual(values = c( "#00AFBB","#FC4E07"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  risk),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)+
  ylim(0,1)+
  coord_flip()
box_TME1; ggsave("TME1.pdf",box_TME1,height=15,width=18,unit="cm")

box_TME2 <- TME_New[TME_New$algorithm == "MCPCOUNTER",] %>%
  mutate(Celltype = factor(Celltype, levels = unique(TME_New$Celltype)))%>% 
  ggplot(aes(x = Celltype, y = log10(Composition), color = risk))+ 
  labs(y="Tumor Environment Score",x= NULL,title = NULL)+  
  geom_boxplot(position=position_dodge(0.5),width=0.5,alpha = 0.8, outlier.size = 0)+ 
  scale_color_manual(values = c( "#00AFBB","#FC4E07"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  risk),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)+
  coord_flip()
box_TME2; ggsave("TME2.pdf",box_TME2,height=6,width=18,unit="cm")


#####################
#TIDE
rm(list=ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
exp<-read.table("TCGA_ov_exp_log2cpm.txt",header=T,sep="\t", check.names = F)
exp <- t(exp)

head(exp[1:4, 1:4])
range(exp)
exp2 <- t(apply(exp, 1, function(x)x-(mean(x)))) 
head(exp2[1:4, 1:4])
range(exp2)

setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/tide/")
write.table(exp2, "tide_exp.txt", col.names = T, row.names = T, sep = "\t")

rm(list=ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/tide/")
tide <- read.csv("tide_total_OV.csv")  

setwd("D:/工作妇肿/课题/disulfidptosis_v2/train_multicox/")
risk <- read.table("coxriskscore.txt", header = T, row.names = 1, check.names = F)
colnames(risk)[1] <- "Patient" 

tide_risk <- merge(tide, risk, by = "Patient")
setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/tide/")
write.csv(tide_risk, "tide_risk.csv",row.names = F, quote = F)

rm(list=ls())
setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/tide/")
tide_risk <- read.csv("tide_risk.csv", check.names = F)

library(ggplot2)
library(ggpubr)
library(ggbeeswarm)

if(T){
  mytheme <-   theme(legend.position="none", plot.title = element_text(hjust = 0.5),
                     axis.text.x = element_text(angle=30, hjust = 1, vjust=1),
                     legend.text = element_text(size= 14),
                     legend.title= element_text(size= 14)
  ) }

var1 <- c("No.benefits", "Responder", "CTL.flag")
var2 <- c("TIDE",  "IFNG", "MSI.Expr.Sig", "Merck18", "CD8", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM.M2")

setwd("D:/工作妇肿/课题/disulfidptosis_v2/immune_analysis/tide/Immune_response/")
outTab=data.frame()
for (i in c("risk")) {
  for(j in var2){
    x <- tide_risk[,i]
    x <- factor(x, levels = c("Low_risk", "High_risk"))
    y <- tide_risk[,j]
    z <- wilcox.test(y~x, paired=F)$p.value
    ptable <- data.frame(i, j, low_median= median(tide_risk[which(tide_risk[,i]=="Low_risk"),j]), 
                         high_median= median(tide_risk[which(tide_risk[,i]=="High_risk"),j]),z)
    outTab <- rbind(ptable, outTab)
    plot <- ggplot(tide_risk, aes(x, y, color = x))+
      geom_violin()+
      geom_jitter(width=0.2, size = 0.9, alpha=0.8)+
      geom_boxplot(stat = "boxplot",
                   fill = "white",
                   width=0.2,
                   alpha = 1,
                   outlier.color = "grey", 
                   outlier.fill = "grey", 
                   outlier.shape = 19,
                   outlier.size = 1, 
                   outlier.alpha = 0.1)+
      stat_compare_means(label = "p.signif", label.x=1.3, hide.ns = F)+ ###ggpubr添加P值和显著性水平
      scale_fill_manual(values = c( "#00AFBB","#FC4E07"))+ 
      scale_color_manual(values = c( "#00AFBB","#FC4E07"))+
      xlab(i)+
      ylab(j)+
      theme_bw()+mytheme
    plot_name <- paste0(i,j,"_boxplot.pdf")
    plot;ggsave(plot_name, plot, height=8,width=5,unit="cm")
  }
}
colnames(outTab)[c(1,2,5)]<- c("gene", "Immune_marker", "p value")
write.csv(outTab, "pvalue_risk.csv")


#scRNA
rm(list = ls())

# 数据来源：PMID 32561858， Sample size 7, cell number 40225
library(Seurat)
library(SingleR)
library(dplyr)
library(patchwork)
library(mindr)
library(Matrix)
library(ggsci)
setwd("E:/细胞死亡模型/scRNA_OV/export/OvC_counts/")
ov.data <- Read10X(data.dir = ".")
ov <- CreateSeuratObject(counts = ov.data, project = "ov3p", min.cells = 3, min.features = 200)
ov

save(ov,file = "ov.rdata") 

rm(list=ls())
load(file = "ov.rdata")

ov[["percent.mt"]] <- PercentageFeatureSet(ov, pattern = "^MT-") #"[["这个操作符可以在列中添加到对象中，使用从MT-作为线粒体基因集
head(ov@meta.data, 5)#在创建对象 CreateSeuratObject() 的过程中，会自动计算细胞中独特基因与总基因数目，可以在目标对象中找到
mean(ov@meta.data$nCount_RNA)

HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
HB_m <- match(HB.genes_total,rownames(ov@assays$RNA))
HB.genes <- rownames(ov@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
ov[["percent.HB"]]<-PercentageFeatureSet(ov,features=HB.genes)

# write  nGene_nUMI_mito_HB
head(ov@meta.data)[,c(2,3,4,5)]

#Feature、count、线粒体基因、红细胞基因占比可视化
pdf("VlnPlots.pdf", height = 10, width = 20)
VlnPlot(ov, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4)  #QC指标可视化 
dev.off()

#均一化和标准化 percent.mt < 10
ov <- subset(ov, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & percent.HB < 10)#从可视化结果进行细胞选择，过滤掉超过6000或少于200基因的细胞，同时过滤线粒体基因比例超过10%、红细胞基因比例超过10%的细胞
table(ov@meta.data$orig.ident)
sum(table(ov@meta.data$orig.ident)) #过滤完后细胞22875


ov <- NormalizeData(ov, normalization.method = "LogNormalize", scale.factor = 10000)#过滤细胞后，进行标准化数据，默认情况下，使用“lognormalize”全局缩放的归一化方法，通过总表达值对每个细胞的基因表达值归一化，并将其乘以缩放因子(默认为10,000)，最后对结果进行对数变换，标准化数据存储在ov@assays$RNA@data中。NormalizeData()这个函数是首先对基因的reads数进行了同一文库大小的校正，然后再对校正后的数值进行对数化操作
ov@assays$RNA@data[1:5,1:5] #NormalizeData()后的数据存储在ov@assays$RNA@data中

#Identification of highly variable features (feature selection)
ov <- FindVariableFeatures(ov, selection.method = "vst", nfeatures = 2000)#计算数据集中表现出高细胞间差异的基因子集(即，它们在一些细胞中高表达，而在另一些细胞中低表达)。在下游分析中关注这些基因有助于突出单细胞数据集中的生物信号。默认情况下，每个数据集选择2000个高变异基因用于下游分析。


# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(ov), 20)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ov)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
pdf("variableFeatures_Top20.pdf", width = 12, height = 10)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
dev.off()

# Scaling the data 
all.genes <- rownames(ov)
ov <- ScaleData(ov, features = all.genes, vars.to.regress = c("percent.mt", "percent.HB") ) #对基因表达量的数值进行了z-score的转换，同时也为后面的pca分析做了铺垫，因为pca分析默认数据是服从正态分布的
ov@assays$RNA@scale.data[1:5, 1:5] #scaled data存储在这里。ScaleData默认使用的基因是FindVariableFeatures找出来的基因
save(ov, file = "ov_0325.rdata")

#内参基因验证是否需要去批次
pdf("gapdh_actb_dotplot.pdf")
DotPlot(ov, features = c("HBB", "ACTB")) 
dev.off()

#Perform linear dimensional reductionPerform linear dimensional reduction
ov <- RunPCA(ov, features = VariableFeatures(object = ov))

head(ov@reductions$pca@cell.embeddings)#每个细胞在PC轴上的坐标
head(ov@reductions$pca@feature.loadings)#每个细胞在PC轴上的坐标


###寻找最佳PC数
# Determine percent of variation associated with each PC
pct <- ov[["pca"]]@stdev / sum( ov[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)
pcs ###算出来最佳pc数为17

# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))

# Elbow plot to visualize 
library(ggplot2)
pdf("bestPCs_elbow.pdf", width = 10, height = 8)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

# Printing out the most variable genes driving PCs
print(x = ov[["pca"]],  dims = 1:17,  nfeatures = 5)

# Cluster the cells 
#Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
ov <- FindNeighbors(ov, dims = 1:17)
ov <- FindClusters(ov, resolution = 0.5)

table(ov@active.ident)

#Calculate the Barcode Distribution Inflection
ov<-CalculateBarcodeInflections(ov)
SubsetByBarcodeInflections(ov)


#Run non-linear dimensional reduction (UMAP/tSNE)
library(uwot)
library(magrittr)
ov <- RunUMAP(ov, dims = 1:17)%>% FindNeighbors(dims = 1:17) %>% FindClusters(resolution = 0.5) %>% identity()
pdf("UMAP_mtDNA10.pdf", width = 10, height = 10)
DimPlot(ov, reduction = "umap",label = T)
dev.off()

p1 <- DimPlot(ov, group.by = "orig.ident")
p2 <- DimPlot(ov,  reduction = "umap", split.by = "orig.ident", ncol = 5)
p2
dev.off()

library(DoubletFinder)
##DoubletFinder
Find_doublet <- function(data){
  #测试最佳参数的过程，运行速度慢
  #使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  #最佳的参数，BCmetric、MeanBC
  #比较重要的一步是找到合适的pK值。
  #作者引入了一个参数BCmetric，是他们通过分析人和老鼠混合单细胞样本（知道真实双细胞率的样本）得出的经验参数，其最大值对应合适的pK值。
  mpK <- as.numeric(as.vector(bcmvn[bcmvn$BCmetric==max(bcmvn$BCmetric),]$pK))
  DoubletRate = ncol(data)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
  nExp_poi <- round(DoubletRate*ncol(data))
  
  #双细胞有两种，同源双细胞和异源双细胞。DoubletFinder只能检测异源双细胞。所以需要把同源双细胞可能的比率去除掉，以优化期望的doublets数量。
  #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞
  #最好提供celltype，而不是seurat_clusters。
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #使用nExp = nExp_poi和nExp = nExp_poi.adj,分别进行doublets鉴定，以便后续确定哪些细胞是Doublet-Low/High Confidience
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_low"
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_high"
  return(data)
}
dim.usage = 17
a <- Find_doublet(ov)

table(ov@meta.data$doublet_low)
table(ov@meta.data$doublet_high)

# Plot results ------------------------
ov@meta.data$DF_hi.lo <- ov@meta.data$doublet_low
ov@meta.data$DF_hi.lo[which(ov@meta.data$DF_hi.lo == "Doublet" & ov@meta.data$doublet_high == "Singlet")] <- "Doublet-Low"
ov@meta.data$DF_hi.lo[which(ov@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High"
table(ov@meta.data$DF_hi.lo)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ov.markers <- FindAllMarkers(ov, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ov.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(ov.markers, "ov.markers.csv")

#DoHeatmap
ov.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10, "ov_top10.csv")

pdf("heatmap_cluster.pdf", width = 25, height = 30)
DoHeatmap(ov, features = top10$gene) + NoLegend()
dev.off()

ov@meta.data$condition <- ov@meta.data$orig.ident
ov@meta.data[,c("condition")] <- gsub("BT1303", "normal", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "BT1304", "normal", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "scrSOL007", "normal", ov@meta.data[,c("condition")])

ov@meta.data[,c("condition")] <- gsub("BT1305", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "BT1306", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "BT1307", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub("scrSOL001", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "scrSOL003", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub( "scrSOL004", "tumor", ov@meta.data[,c("condition")])
ov@meta.data[,c("condition")] <- gsub("scrSOL006", "tumor", ov@meta.data[,c("condition")])
table(ov@meta.data$condition)

##cell annotation
pdf("dotplot_umap_marker_seurat_clusters_mtdna10.pdf", width = 10, height = 10)
DotPlot(ov, 
        assay = ov@reductions$pca@assay.used,
        features = c("LYZ", "CD68", "CD14", "CD86", "HLA-DRA", "FCGR3A", #myeloid
                     "EPCAM", "KRT7", "KRT8", "KRT17", "SPRR3", #tumor/epithelial
                     "CD3E", "CD3D", "TRBC1", "TRBC2", "TRAC", #t cell
                     "CD79A", "CD79B", "JCHAIN", "IGKC", "IGHG3", #B cell/plasma cell
                     "CLDN5", "FLT1", "CDH1", "RAMP2", #endothelial cell
                     "DCN", "C1R", "COL1A1", "ACTA2", #fibroblast
                     "TAGLN", "CNN1", #Smooth muscle cell
                     "TPSAB1" #mast cell
        )
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

##cell annotation reference：https://doi.org/10.1038/s41467-021-27599-5

celltype <- data.frame(clusterID = 0:22,
                       celltype = "other")

celltype[celltype$clusterID %in% c(2, 9, 16, 21), 2] <- "myeloid cell"
celltype[celltype$clusterID %in% c(6, 7, 10, 11, 14), 2] <- "tumor/epithelial cell"
celltype[celltype$clusterID %in% c(1, 19), 2] <- "T cell"
celltype[celltype$clusterID %in% c(15, 20), 2] <- "B cell/plasma cell"
celltype[celltype$clusterID %in% c(8), 2] <- "endothelial"
celltype[celltype$clusterID %in% c(0, 3, 4, 5, 12, 13, 18, 22), 2] <- "fibroblast"
celltype[celltype$clusterID %in% c(17), 2] <- "smooth muscle cell"

table(celltype$celltype)

for( i in 1:nrow(celltype)){
  ov@meta.data[which (ov@meta.data$seurat_clusters == celltype$clusterID[i]), "celltype"] <- 
    celltype$celltype[i]
}
table(ov@meta.data$celltype)
write.csv(ov@meta.data, "ov1_metadata.csv")
saveRDS(ov,file = "ov_celltype_0325.RDS") #.rdata

ov <- SetIdent(ov, value = "celltype")
ov@active.ident <- factor(ov@active.ident, 
                              levels=c("myeloid cell", 
                                       "tumor/epithelial cell",
                                       "T cell",
                                       "B cell/plasma cell",
                                       "endothelial",
                                       "fibroblast",
                                       "smooth muscle cell"
                              ))

pdf("dotplot_umap_marker_celltype_0325.pdf", width = 10, height = 3)
DotPlot(ov, 
        assay = ov@reductions$umap@assay.used,
        features = c("LYZ", "CD68", "CD14", "CD86", "HLA-DRA", "FCGR3A", #myeloid
                     "EPCAM", "KRT7", "KRT8", "KRT17", "SPRR3", #tumor/epithelial
                     "CD3E", "CD3D", "TRBC1", "TRBC2", "TRAC", #t cell
                     "CD79A", "CD79B", "JCHAIN", "IGKC", "IGHG3", #B cell/plasma cell
                     "CLDN5", "FLT1", "CDH1", "RAMP2", #endothelial cell
                     "DCN", "C1R", "COL1A1", "ACTA2", #fibroblast
                     "TAGLN", "CNN1", #Smooth muscle cell
                     "TPSAB1" #mast cell
        )
)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

ov@meta.data$seurat_clusters
library(ggsci)
library(cowplot)
p1 <- DimPlot(ov, reduction = "umap",  group.by = "celltype", label = TRUE)+ scale_color_npg()
p8 <- DimPlot(ov, reduction = "umap",  group.by = "seurat_clusters", label = TRUE)
p2 <- FeaturePlot(ov, "MYL6", min.cutoff = 1)
p3 <- FeaturePlot(ov, "PDLIM1")
p4 <- FeaturePlot(ov, "ACTN4")
p5 <- FeaturePlot(ov, "FLNB")
p6 <- FeaturePlot(ov, "SLC7A11")
p7 <- FeaturePlot(ov, "CD2AP")

pdf("disulfidpsis_marker_0325.pdf", width = 10, height = 4)
p9 <- p8 + p1
plot_grid(p9)
dev.off()

pdf("disulfidpsis_markerGenes_0325.pdf", width = 10, height = 6)
p10 <- p2 + p3 + p4 + p5 + p6 + p7
plot_grid(p2, p3, p4, p5, p6, p7, ncol = 3)
dev.off()

library(ggplot2)
p1 <- VlnPlot(ov, features = c("MYL6"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))
p2 <- VlnPlot(ov, features = c("PDLIM1"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))
p3 <- VlnPlot(ov, features = c("ACTN4"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))
p4 <- VlnPlot(ov, features = c("FLNB"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))
p5 <- VlnPlot(ov, features = c("SLC7A11"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))
p6 <- VlnPlot(ov, features = c("CD2AP"),group.by = "celltype", pt.size = 0, split.by = "condition" ) + 
  scale_fill_manual(values= c( "tumor" = "orange", "normal" = "blue"))

pdf("disulfidpsis_markerGenes_tumor.vs.normal.0325.pdf", width = 15, height = 10)
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
dev.off()

##expression distribution of MYL6
condition <- ov@meta.data$condition
celltype <- ov@meta.data$celltype
MYL6 <- as.numeric(ov@assays$RNA["MYL6",])

df <- data.frame(MYL6=MYL6,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("MYL6_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "MYL6", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()

##expression distribution of PDLIM1
PDLIM1 <- as.numeric(ov@assays$RNA["PDLIM1",])

df <- data.frame(PDLIM1=PDLIM1,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("PDLIM1_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "PDLIM1", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()

##expression distribution of ACTN4
ACTN4 <- as.numeric(ov@assays$RNA["ACTN4",])

df <- data.frame(ACTN4=ACTN4,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("ACTN4_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "ACTN4", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()

##expression distribution of FLNB
FLNB <- as.numeric(ov@assays$RNA["FLNB",])

df <- data.frame(FLNB=FLNB,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("FLNB_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "FLNB", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()

##expression distribution of SLC7A11
SLC7A11 <- as.numeric(ov@assays$RNA["SLC7A11",])

df <- data.frame(SLC7A11=SLC7A11,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("SLC7A11_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "SLC7A11", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()


##expression distribution of CD2AP
CD2AP <- as.numeric(ov@assays$RNA["CD2AP",])

df <- data.frame(CD2AP=CD2AP,
                 condition=condition,
                 celltype = celltype)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12))

library(dplyr)

## `summarise()` ungrouping output (override with `.groups` argument)
df$celltype <-  factor(df$celltype)

library(ggpubr)
library(ggplot2)

pdf("CD2AP_tumor_vs._normal_myeloid.pdf",height=5,width=8)
ggboxplot(df, "celltype", "CD2AP", color = "condition", fill= "white", group= "celltype", )+
  scale_color_manual(values = c( "#2546F1", "#E4392E"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  condition),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
dev.off()


# TCGA vs. GTEx
rm(list = ls())

setwd("./TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2/")

library(data.table)
gtex_exp<-fread("TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2", header = T,sep = "\t", check.names = F)
head(gtex_exp[1:4,1:4])
row.names(gtex_exp) <- gtex_exp$V1

pheno <- fread("TcgaTargetGTEX_phenotype.txt", header = T,sep = "\t", check.names = F)
head(pheno[1:4,1:4])
row.names(pheno) <- pheno$Sample
a <- unique(pheno$'_primary_site')
length(a)
print(a)

library(magrittr)
library(tidyverse)
library(dplyr)
b <- pheno[which(pheno$'_primary_site'=="Ovary"),]$sample 
c <- colnames(gtex_exp)
d <- intersect(b,c)
d<- append("V1", d)
ov <- gtex_exp %>% select(d)
head(ov[1:4,1:4])

ano<-fread("probeMap_gencode.v23.annotation.gene.probemap", header = T,sep = "\t", check.names = F)
colnames(ov)[1]<- "id"
ov_ano <- merge(ano[,1:2],ov, by="id")
head(ov_ano[1:4,1:10])

ov_ano <- ov_ano[,-1]

library(limma)
ov_ano_2=avereps(ov_ano,ID=ov_ano$gene) #对重复的基因的表达量取均值
dim(ov_ano_2)

write.table(ov_ano_2,"ov_ano_2.txt",row.names = F,col.names = T,sep = "\t")

rm(list = ls())
########load GTEx data
library(data.table)
ov_ano<-fread("ov_ano_2.txt", header = T,sep = "\t", check.names = F)
head(ov_ano[1:4,1:4])
gene<-ov_ano$gene
row.names(ov_ano) <- gene
ov_ano <- ov_ano[,-1]
dim(ov_ano)
head(ov_ano[1:4,1:4])
ov_ano<-sapply(ov_ano,as.numeric,2)
ov_ano<- 2^(ov_ano)-1 ##将GTEx数据从log2(count+1)转化为count格式
rownames(ov_ano)<- gene
head(ov_ano[1:4,1:4])  ##ov_ano是gtex count形式的表达矩阵
ov_ano <- round(ov_ano,0)

ov_tcga <- ov_ano[, grep("TCGA", colnames(ov_ano))]  ###ov_ano里TCGA的数据
dim(ov_tcga)

c4<-sapply(strsplit(colnames(ov_tcga),"-"),"[",4)
c5<-sapply(strsplit(c4,""),"[",1)
table(c5) ##########发现TCGA的427个人全是0开头的（全是tumor）。那就以TCGA和GTEx分tumor和normal

#然后因为这个表是Xena已经合并且标准化过的，所以也不用再去批次啦
write.table(ov_ano, "ov_TCGA_GTEx_count.txt",row.names = T,col.names = T,sep = "\t")

rm(list = ls())
ov_ano<-read.table("ov_TCGA_GTEx_count.txt", header = T,sep = "\t", check.names = F)
multiInput=read.csv("multiCox.csv")

#limma包分析差异表达基因
library(limma)
ov_pred <- ov_ano[which(rownames(ov_ano) %in% multiInput$id),]  #TCGA_GTEx 纳入的28个基因
ov_pred <- log2(ov_pred+0.1)
list <- model.matrix(~factor(condition)+0)
table(list)
colnames(list) <- c("GTEX", "TCGA")
df.fit <- lmFit(ov_pred,list)

df.matrix <- makeContrasts(TCGA-GTEX, levels=list)
fit<- contrasts.fit(df.fit,df.matrix)
fit<- eBayes(fit)
tempOutput<- topTable(fit, n=Inf, adjust= "fdr")
head(tempOutput)

nrDEG = na.omit(tempOutput) 
write.csv(diffsig, "all.limmaOut.csv")

foldChange = 2
padj = 0.001
All_diffSig <- diffsig[(diffsig$padj < padj & (diffsig$log2FoldChange>foldChange | diffsig$log2FoldChange < (-foldChange))),]
dim(All_diffSig)
write.csv(All_diffSig, "all.diffsig_tumor_normal.csv")  

diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$log2FoldChange > foldChange)),]
write.csv(diffup, "diffup.csv")

diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$log2FoldChange < -foldChange)),]
write.csv(diffdown, "diffdown.csv")

###Tumor vs. normal valcano plot
library(ggplot2)
library(ggrepel)
logFC <- diffsig$log2FoldChange
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$log2FoldChange < foldChange) & data$log2FoldChange > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$log2FoldChange > 1)] <-  "Up"
data$group[(data$padj <= 0.05 & data$log2FoldChange < -1)] <- "Down"
x_lim <- max(logFC,-logFC)

pdf('tumor_normal_volcano.pdf',width = 7,height = 6.5) 
label = subset(diffsig,P.Value <0.05 & abs(logFC) > 0.5)
label1 = rownames(label)

colnames(diffsig)[1] = 'log2FC'
Significant=ifelse((diffsig$P.Value < 0.05 & abs(diffsig$log2FC)> 0.5), ifelse(diffsig$log2FC > 0.5,"Up","Down"), "Not")

ggplot(diffsig, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig, max.level = c(-1, 1))+theme_bw()

dev.off()

#DRPS functional enrichment analysis
library(clusterProfiler)

gene_list = c("MYL6", "PDLIM1", "ACTN4", "FLNB", "SLC7A11", "CD2AP")

gene.df <- bitr(gene_list, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
gene.df <- na.omit(gene.df)

GO <- enrichGO(gene = gene.df$ENTREZID, 
               OrgDb = "org.Hs.eg.db", 
               pvalueCutoff = 0.05,
               ont="all")

a <- GO@result

library(ggplot2)
setwd("D:/工作妇肿/课题/disulfidptosis_v2/GSEA/")
pdf("GO_dotplot_drps.pdf", width = 8, height = 16)
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dev.off()

#KEGG
KEGG_database <- 'hsa' 
KEGG<-enrichKEGG(gene.df$ENTREZID,
                 organism = KEGG_database,
                 pvalueCutoff = 0.5,
                 minGSSize = 10, maxGSSize = 5000
)

library(org.Hs.eg.db)
ekegg <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_diff_dt <- data.frame(ekegg)
p1 <- dotplot(ekegg, showCategory=20, title = 'KEGG Pathway')

pdf("KEGG_dotplot_drps.pdf", width = 8, height = 6)
p1
dev.off()


###correlation of DRPS
rm(list = ls())
multiInput_2 <- read.table("risk_train.txt",header=T,sep="\t")
exp <- multiInput_2[,4:9]

library(Hmisc)
res2 <- rcorr(as.matrix(exp), type="pearson")

library(corrplot)
pdf("drps_correlation.pdf")
corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.05,insig = "blank")
dev.off()


#sensitivity analysis
library(compound.Cox)
rm(list = ls())
multiInput_2 <- read.table("risk_train.txt",header=T,sep="\t")
CG.test(multiInput_2$fustime, (multiInput_2$fustat-1), multiInput_2$riskScore, 
        copula= CG.Frank,alpha=17, N=60)

multiInput_2 <- read.table("risk.GSE9891.txt",header=T,sep="\t")
multiInput_2 <- na.omit(multiInput_2)
CG.test(multiInput_2$fustime, (multiInput_2$fustat), multiInput_2$riskScore, 
        copula= CG.Clayton,alpha=10, N=60)

multiInput_2 <- read.table("risk.E_MTAB_386.txt",header=T,sep="\t")
multiInput_2 <- na.omit(multiInput_2)
CG.test(multiInput_2$fustime, (multiInput_2$fustat), multiInput_2$riskScore, 
        copula= CG.Frank ,alpha=17, N=60)

#IC50 prediction
rm(list = ls())  
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
setwd("G:/ov/disulfidptosis/DataFiles/")
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 


testExpr<- read.table("G:/ov/disulfidptosis/TCGA_ov_exp_log2cpm.txt", sep = "\t", check.names = F)
testExpr[1:4,1:4]  
testExpr <- t(testExpr)
testExpr <- testExpr[which(row.names(testExpr) %in% row.names(GDSC2_Expr)), ]
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]
row.names(testPtype) <- testPtype$V1
testPtype <- testPtype[,-1]
write.csv(testPtype, "IC50_ov.csv")

testPtype<- cbind(id=rownames(testPtype), testPtype)

setwd("G:/ov/disulfidptosis/")
risk <- read.table("coxriskscore.txt")
ic50_risk <- cbind(risk, testPtype, by= "id")

df <- data.frame()
drugs <- colnames(testPtype)[c(7:ncol(testPtype))]
for (i in drugs){
  df[i,1] = t.test(as.numeric(ic50_risk[,i])~
                     ic50_risk$risk)$p.value
  df[i,2] = median(ic50_risk[which(ic50_risk$risk=="High_risk"),i])
  df[i,3] = median(ic50_risk[which(ic50_risk$risk=="Low_risk"),i])
}

colnames(df) <- c("pvalue", "median_highRisk", "median_lowRisk")
df <- cbind(drug=rownames(df), df)
df_sig <-df[df$pvalue<0.05 & df$median_highRisk>df$median_lowRisk,]
write.csv(df_sig, "lowRisk_sensitive.csv")

library(ggplot2)
library(ggsignif)

df_sig2 <-df[df$pvalue<0.05 & df$median_highRisk<df$median_lowRisk,]
write.csv(df_sig2, "highRisk_sensitive.csv")

a <- as.data.frame(cbind(risk= ic50_risk$risk, drug = ic50_risk$`Staurosporine_1034`))
a$drug <- as.numeric(a$drug)
pdf("Staurosporine_1034.pdf", width = 4, height = 4)
ggplot(a,aes(x=risk, y=drug, color = risk))+
  geom_boxplot(width=0.5)+
  #geom_jitter(shape=16, position = position_jitter(0.2))
  theme_bw()+
  geom_signif(comparisons = list(c("High_risk", "Low_risk")),
              map_signif_level=TRUE)
dev.off()

a <- as.data.frame(cbind(risk= ic50_risk$risk, drug = ic50_risk$`Tamoxifen_1199`))
a$drug <- as.numeric(a$drug)
pdf("Tamoxifen_1199.pdf", width = 4, height = 4)
ggplot(a,aes(x=risk, y=drug, color = risk))+
  geom_boxplot(width=0.5)+
  #geom_jitter(shape=16, position = position_jitter(0.2))
  theme_bw()+
  geom_signif(comparisons = list(c("High_risk", "Low_risk")),
              map_signif_level=TRUE)
dev.off()

a <- as.data.frame(cbind(risk= ic50_risk$risk, drug = ic50_risk$`Epirubicin_1511`))
a$drug <- as.numeric(a$drug)
pdf("Epirubicin_1511.pdf", width = 4, height = 4)
ggplot(a,aes(x=risk, y=drug, color = risk))+
  geom_boxplot(width=0.5)+
  #geom_jitter(shape=16, position = position_jitter(0.2))
  theme_bw()+
  geom_signif(comparisons = list(c("High_risk", "Low_risk")),
              map_signif_level=TRUE)
dev.off()

a <- as.data.frame(cbind(risk= ic50_risk$risk, drug = ic50_risk$`Navitoclax_1011`))
a$drug <- as.numeric(a$drug)
pdf("Navitoclax_1011.pdf", width = 4, height = 4)
ggplot(a,aes(x=risk, y=drug, color = risk))+
  geom_boxplot(width=0.5)+
  #geom_jitter(shape=16, position = position_jitter(0.2))
  theme_bw()+
  geom_signif(comparisons = list(c("High_risk", "Low_risk")),
              map_signif_level=TRUE)
dev.off()



#############################################################
## GSE18920 --- Splicing array data from human spinal cord ##
#############################################################
library(limma)
path = "/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/"
setwd(path)
data_spinal = read.csv(sprintf("%sGSE18920_extended_for_R_wGeneNames.txt",path),stringsAsFactors = FALSE,sep="\t")
rownames(data_spinal) = gsub(" ","",data_spinal$NAME)
data_spinal = data_spinal[,-1]
data_spinal = data_spinal[,grep("AH",colnames(data_spinal))] #<--- THIS IS THE LINE I SHOULD HAVE HAD PREVIOUSLY
annot_spinal_status = rep("ALS",dim(data_spinal)[2])
annot_spinal_status[grep("CTRL",colnames(data_spinal))] = "CTRL"
annot_spinal_gender = rep("Male",dim(data_spinal)[2])
annot_spinal_gender[grep("_F",colnames(data_spinal))] = "Female"
annot_spinal_age = as.numeric(as.character(gsub("M|F","",gsub("No.*","",gsub(".*_","",colnames(data_spinal))))))
annot_spinal = data.frame(status=annot_spinal_status,gender=annot_spinal_gender,age=annot_spinal_age)
annot_spinal_tissue = gsub("_.*","",colnames(data_spinal))

prep.tt <- function(exp_data,mod,m2h,title,coef){
    fit = lmFit(exp_data,mod)
    eb = eBayes(fit)
    tt = topTable(eb, coef=coef, adjust="BH",number=1000000)	
    tt2 = cbind(tt,HGNC.symbol=rownames(tt))
    tt2$HGNC.symbol=as.character(tt2$HGNC.symbol)
    tt3 = merge(tt2,m2h,by="HGNC.symbol")
    tt4 = tt3[order(tt3$P.Value),]
    #colnames(tt4)[length(colnames(tt4))]="HGNC.symbol"
    #tt_mus = merge(tt4,m2h,by="HGNC.symbol")
    #tt_mus2 = tt_mus[order(tt_mus$P.Value),]
    #print(sprintf("TopTables/%s",title))
    #colnames(tt_mus2)[9]="Gene_symbol"
    write.csv(tt4,file=sprintf("TopTables/%s",title))
    return(tt4)
}

load("/Users/natske/Google Drive/DiseaseEnrichment/mouse_to_human_homologs.rda")
m2h = unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])

mod  = model.matrix(~annot_spinal$gender+factor(annot_spinal$status,levels=c("CTRL","ALS")))
colnames(mod)[2:3] = c("Gender","ALS")
#mod  = model.matrix(~annot_spinal$age+annot_spinal$gender+factor(annot_spinal$status,levels=c("CTRL","ALS")))
#colnames(mod)[2:4] = c("Age","Gender","ALS")
#mod  = model.matrix(~factor(annot_spinal$status,levels=c("CTRL","ALS"))+annot_spinal_tissue+annot_spinal_tissue)
tt_spinal = prep.tt(data_spinal,mod,m2h,"tt_HumanSpinal.csv",coef="ALS")

### NOW RUN EWCE

source("/Users/natske/Google Drive/EWCE_old/R/generate.bootstrap.plots.r")
source("/Users/natske/Google Drive/EWCE_old/R/read_celltype_data.r")
source("/Users/natske/Google Drive/EWCE_old/R/get_annot_below.r")
source("/Users/natske/Google Drive/EWCE_old/R/bootstrap.enrichment.test.r")
source("/Users/natske/Google Drive/EWCE_old/R/ewce.plot.r")
source("/Users/natske/Google Drive/EWCE_old/R/merge_two_expfiles.r")
source("/Users/natske/Google Drive/EWCE_old/R/ewce_expression_data.r")
source("/Users/natske/Google Drive/EWCE_old/R/get_summed_proportions.r")
source("/Users/natske/Google Drive/EWCE_old/R/cell.list.dist.r")

thresh=0#.00
trim=0.0
wtLEVELS="woLev"
sct_path = "/Users/natske/Datasets that are too large to store elsewhere/SOD1 Spinal Cord (GSE18597)"
#load(file=sprintf("%s/celltype_data_OligosNCortex_(%s)_thresh(%s)_trim(%s).rda",sct_path,wtLEVELS,thresh,trim))
load("/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/celltype_data_OligosNCortex_(woLev)_thresh(0)_trim(0).rda")

#thresh = 100
thresh = 250
#sct_data = celltype_data
#annot=celltype_data[[1]]$annot
#annotLevel=1
#tt=tt_spinal
#sortBy="t"
#reps=100
#useHGNC=FALSE

#save(tt_spinal,file="tt_spinal.rda")
#load(file="tt_spinal.rda")

ewce_spinal = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_spinal,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
ewce_spinal$joint_results$Q = p.adjust(ewce_spinal$joint_results$p)
ewce_spinal$joint_results[order(ewce_spinal_NEW$joint_results$sd_from_mean),]
#save(ewce_spinal_NEW,file="EWCE_SPINAL.Rda")
load(ewce_spinal,file="EWCE_SPINAL.Rda")
write.csv(ewce_spinal$joint_results,file="EWCE_SPINAL.csv")
pdf(file="Fig_EWCE_HumanSpinal.pdf",width=6,height=5)
ewce.plot(ewce_spinal$joint_results)
dev.off()

source("/Users/natske/Google Drive/EWCE/R/generate.bootstrap.plots.for.transcriptome.r")
#load("/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/EWCE_SPINAL.Rda")
generate.bootstrap.plots.for.transcriptome(sct_data=celltype_data,tt=tt_spinal,thresh=100,annotLevel=1,reps=1000,full_results=ewce_spinal,listFileName="Human ALS Spinal")

#############################################################
## GSE18920 --- Splicing array data from human spinal cord ##
#############################################################
library(limma)
path = "/Users/natske/Datasets that are too large to store elsewhere/Disease Transcriptomes/ALS/"
setwd(path)
exp_cortex = read.csv(sprintf("%sGSE67196_Petrucelli2015_ALS_genes.rawcount_EXP.csv",path),stringsAsFactors = FALSE,sep=",")
exp_cortex = exp_cortex[!duplicated(exp_cortex$GeneID),] 
rownames(exp_cortex) = exp_cortex$GeneID
exp_cortex = exp_cortex[,-c(1:2)]

annot_cortex =  read.csv(sprintf("%sGSE67196_Petrucelli2015_ALS_genes.rawcount_ANNOT.csv",path),stringsAsFactors = FALSE,sep=",")
#annot_cortex$X.Clinical.diagnosis[annot_cortex$X.Clinical.diagnosis=="aMCI"]="Healthy"
annot_cortex$X.Clinical.diagnosis[annot_cortex$X.Clinical.diagnosis=="Healthy (depression)"]="Healthy"
annot_cortex$X.Clinical.diagnosis[annot_cortex$X.Clinical.diagnosis=="Bulbar palsy"]="ALS"
annot_cortex$X.Clinical.diagnosis[annot_cortex$X.Clinical.diagnosis=="ALS-MCI"]="ALS"
annot_cortex        = annot_cortex[annot_cortex$X.Clinical.diagnosis %in% c("Healthy","ALS"),]
annot_cortex_status = annot_cortex$X.Clinical.diagnosis
annot_cortex_age    = annot_cortex$Age.of.death.years.
annot_cortex_gender = annot_cortex$sex
annot_cortex        = data.frame(group=annot_cortex$Group,sampleNo=annot_cortex$Sample..,status=annot_cortex_status,gender=annot_cortex_gender,age=annot_cortex_age)

sample_names = rep("",dim(annot_cortex)[1])
sample_names[as.character(annot_cortex$group)=="C9ALS"]="C9ALS"
sample_names[as.character(annot_cortex$group)=="sALS"]="sALS"
sample_names[as.character(annot_cortex$group)=="Non-neurological disease control"]="ContrALS"
library(stringr)
sample_names = sprintf("%s%s",sample_names,str_pad(annot_cortex$sampleNo, 3, pad = "0"))
sample_names_fcx = sprintf("%s_fcx",sample_names)
sample_names_cereb = sprintf("%s_cereb",sample_names)
annot_cortex$fcx = sample_names_fcx
annot_cortex$cereb = sprintf("ALS%s_cereb",str_pad(annot_cortex$sampleNo, 3, pad = "0"))
annot_fcx = annot_cortex[annot_cortex$fcx %in% colnames(exp_cortex),]
exp_fcx = exp_cortex[,annot_fcx$fcx]
exp_fcx = exp_fcx[apply(exp_fcx,1,sum)>0,]
annot_cereb = annot_cortex[annot_cortex$cereb %in% colnames(exp_cortex),]
exp_cereb = exp_cortex[,annot_cereb$cereb]
exp_cereb = exp_cereb[apply(exp_cereb,1,sum)>0,]

mod  = model.matrix(~annot_fcx$age+annot_fcx$gender+factor(annot_fcx$status,levels=c("Healthy","ALS")))
colnames(mod)[2:4] = c("Age","Gender","ALS")
tt_fcx = prep.tt(exp_fcx,mod,m2h,"tt_HumanFrontalCortex.csv",coef="ALS")


mod  = model.matrix(~annot_cereb$age+annot_cereb$gender+factor(annot_cereb$status,levels=c("Healthy","ALS")))
colnames(mod)[2:4] = c("Age","Gender","ALS")
tt_cereb = prep.tt(exp_cereb,mod,m2h,"tt_HumanCereb.csv",coef="ALS")

ewce_fcx = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_fcx,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
save(ewce_fcx,file="EWCE_FCX.Rda")
write.csv(ewce_fcx$joint_results,file="EWCE_FCX.csv")
pdf(file="Fig_EWCE_HumanFrontalCortex.pdf",width=6,height=5)
ewce.plot(ewce_fcx$joint_results)
dev.off()


ewce_cereb = ewce_expression_data(celltype_data,annot=celltype_data[[1]]$annot,annotLevel=1,tt_cereb,sortBy="t",thresh=thresh,reps=10000,useHGNC=FALSE)
save(ewce_cereb,file="EWCE_CEREB.Rda")
write.csv(ewce_cereb$joint_results,file="EWCE_CEREB.csv")
pdf(file="Fig_EWCE_HumanCereb.pdf",width=6,height=5)
ewce.plot(ewce_cereb$joint_results)
dev.off()
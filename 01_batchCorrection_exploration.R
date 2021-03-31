# projects to consider
tumors_considered <- c("BLCA","BRCA","CHOL","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PRAD","THCA","UCEC")


# Check NAs in callrates



# Check callRate pval
 pvals <- dpval(result$rnb.set, row.names=TRUE)
x1=apply(pvals,2,is.na)
x2=colSums(x1)


# explore batch correction for 3 tumor types: BRCA, KIRP and 

#
system("mkdir /root/TCGA/ex3")
system("mkdir /root/TCGA/ex3/Rnbeads")
setwd("/root/TCGA/ex3/Rnbeads")
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
suppressMessages(library(IlluminaHumanMethylation450kmanifest))

## preprocessing
#idat files
idat.dir <- file.path("/root/TCGA/tcgaBiolink/")

#### Link list processing ####
# master link list
link.list <- read.table("/root/TCGA/tcgaBiolink/idat_filename_case.txt",header=T,sep="\t")
# remove double associated methylation profiles
link.list <- link.list[grep(",",link.list$cases,invert=TRUE),]
# get project names
projects <- sort(unique(gsub("TCGA\\-","",link.list$project,perl=TRUE)))
# Make sure it matches the panCan data
panCan.dir <- list.dirs(path = "/root/TCGA/panCancer_2018", full.names = TRUE, recursive = FALSE)
for(i in 1:length(projects)){  if(grep(tolower(projects[i]),panCan.dir) > 0){print(paste(projects[i],": OK"))}  }

# Get all genes that can have mutations across all projects
master.gene.list <- read.table(pipe("cat /root/TCGA/panCancer_2018/*/data_mutations_mskcc.txt|tail -n +2|cut -f1|sort|uniq"),sep="\n")
master.gene.list <- as.character(master.gene.list[,1])
##############################

i = 13 #KIRP

  #project specific link.list
  print(paste("Parsing",projects[i],"Files"))
  
  i.link.list = link.list[grep(projects[i],link.list$project),]
  

  #### subset samples with mutation information ####
  mut.file <- paste(panCan.dir[grep(tolower(projects[i]),panCan.dir)],"/data_mutations_mskcc.txt",sep="")
  mut.file <- paste("cut -f1,10,17,40 ",mut.file)
  mut.file <- read.table(pipe(mut.file),sep="\t",header=T,quote="")
  
  meth.id.original <- unique(as.character(i.link.list$cases))
  meth.id <- data.frame( do.call( rbind, strsplit( meth.id.original, '-' ) ) )
  sample.type <- gsub("\\w$","",meth.id[,4],perl=TRUE)
  patient.id <- paste(meth.id[,1],meth.id[,2],meth.id[,3],sep="-")
  meth.id <- paste(meth.id[,1],meth.id[,2],meth.id[,3],sample.type,sep="-")
  
  meth.id.withMutation <- meth.id.original[meth.id %in% mut.file$Tumor_Sample_Barcode]
  mut.id.withMutation <- meth.id[meth.id %in% mut.file$Tumor_Sample_Barcode]
  patient.id.withMutation <- patient.id[meth.id %in% mut.file$Tumor_Sample_Barcode]
  meth.id.normal <- meth.id.original[as.numeric(sample.type)>9 & as.numeric(sample.type)<20]
  
  # create project directory
  system(paste("mkdir",projects[i]))
  
  # Original mutation matrix
  mut.file.ix <- mut.file[mut.file$Tumor_Sample_Barcode %in% mut.id.withMutation,]
  write.csv(mut.file.ix,paste(projects[i],"/",projects[i],"_mutation_original_table.csv",sep=""),row.names=F) #writeStep
  
  # Mutation Matrix
  mut.file.ix <- mut.file.ix[mut.file.ix$Variant_Classification %in% c("Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation",
                                            "Nonstop_Mutation","In_Frame_Del","Frame_Shift_Ins","Nonstop_Mutation"),]
  
  mut.mat <- matrix(0, nrow=length(meth.id.withMutation), ncol=length(master.gene.list) )
  rownames(mut.mat) <- meth.id.withMutation
  colnames(mut.mat) <- master.gene.list
  mut.mat.ix <- unique(mut.file.ix[,c(1,3)])  
  for(ix in 1:dim(mut.mat.ix)[1]){ mut.mat[ mut.id.withMutation %in% mut.mat.ix[ix,2] , master.gene.list %in% mut.mat.ix[ix,1] ] = 1 }
  write.csv(mut.mat,paste(projects[i],"/",projects[i],"_mutation_matrix.csv",sep=""),row.names=F) #writeStep
  if(length(meth.id.normal)>0){ tmp <- matrix(5, nrow=length(meth.id.normal), ncol=length(master.gene.list) )
                                rownames(tmp) <- colnames(meth.id.normal); 
                                colnames(tmp) <- colnames(master.gene.list); 
                               mut.mat <- rbind( mut.mat, tmp) }
  write.csv(mut.mat,paste(projects[i],"/",projects[i],"_mutation_matrix_withNormal.csv",sep=""),row.names=F) #writeStep
  
  # TP53 information matrix
  mut.file.p53 <- mut.file.ix[mut.file.ix$Hugo_Symbol=="TP53",]
  mut.file.p53 <- unique(mut.file.p53)
  
  mut.file.p53 <- cbind(aggregate(as.character(Variant_Classification) ~ 
                     as.character(Tumor_Sample_Barcode), data = mut.file.p53,paste, collapse = "||"),
           aggregate(as.character(HGVSp_Short) ~ 
            as.character(Tumor_Sample_Barcode), data = mut.file.p53,paste, collapse = "||") )
  
  mut.file.p53 <- data.frame(Tumor_Sample_Barcode=mut.file.p53[,1], Variant_Classification=mut.file.p53[,2],HGVSp_Short=mut.file.p53[,4] )
  mut.file.p53 <- data.frame(meth.id.withMutation,mut.file.p53[match(mut.id.withMutation, as.character(mut.file.p53[,1]) ),])
  mut.file.p53 <- mut.file.p53[,c(1,3,4)]
  mut.file.p53[,2] <- as.character(mut.file.p53[,2])
  mut.file.p53[,2][is.na(mut.file.p53[,2])] <- "WT"
  mut.file.p53[,3] <- as.character(mut.file.p53[,3])
  mut.file.p53[,3][is.na(mut.file.p53[,3])] <- "WT"
  write.csv(mut.file.p53,paste(projects[i],"/",projects[i],"_TP53_mutation_info.csv",sep=""),row.names=F) #writeStep
  if(length(meth.id.normal)>0){ tmp <- cbind(meth.id.normal,"NORMAL","NORMAL");
                                colnames(tmp) <- colnames(mut.file.p53); 
                               mut.file.p53 <- rbind( mut.file.p53, tmp) }
  rownames(mut.file.p53) <- NULL
  write.csv(mut.file.p53,paste(projects[i],"/",projects[i],"_TP53_mutation_info_withNormal.csv",sep=""),row.names=F) #writeStep
  
  # CNA table
  CNA <- paste(panCan.dir[grep(tolower(projects[i]),panCan.dir)],"/data_CNA.txt",sep="")
  #CNA <- paste(panCan.dir[i],"/data_CNA.txt",sep="")
  CNA <- read.table(CNA,sep="\t",header=TRUE,quote="")
  CNA <- CNA[!duplicated(CNA$Hugo_Symbol), ]
  rownames(CNA) <- CNA$Hugo_Symbol
  CNA <- CNA[,3:dim(CNA)[2]]
  CNA <- t(CNA)
  rownames(CNA) <- gsub("\\.","\\-",rownames(CNA),perl=TRUE)
  CNA <- CNA[rownames(CNA) %in% mut.id.withMutation,]
  write.csv(CNA,paste(projects[i],"/",projects[i],"_CNA_matrix.csv",sep="")) #writeStep

  # Clinical table
  clinical <- paste(panCan.dir[grep(tolower(projects[i]),panCan.dir)],"/data_clinical_patient.txt",sep="")
  #clinical <- paste(panCan.dir[i],"/data_clinical_patient.txt",sep="")
  clinical <- paste("tail -n+5",clinical,"|cut -f 1,5,6,7,8,9,10,11,26,27,28,29,30")
  clinical <- read.table(pipe(clinical),sep="\t",header=T,quote="")
  clinical.patient.id <- data.frame( do.call( rbind, strsplit( mut.id.withMutation, '-' ) ) )
  clinical.patient.id <- paste(clinical.patient.id[,1],clinical.patient.id[,2],clinical.patient.id[,3],sep="-")
  clinical <- clinical[match(clinical.patient.id, as.character(clinical$PATIENT_ID) ),]
  clinical$PATIENT_ID <- meth.id.withMutation
  write.csv(clinical,paste(projects[i],"/",projects[i],"_clinical_info.csv",sep=""),row.names=F) #writeStep

  ##############################
  print(paste("Starting 450K processing of",projects[i]))


  # Sample annotation creation
  master.meth.id <- c(meth.id.withMutation,meth.id.normal)
  sentrix <- as.character(i.link.list$file_name)
  sentrix <- gsub("_Red.idat","",sentrix)
  sentrix <- gsub("_Grn.idat","",sentrix)
  sentrix <- data.frame( do.call( rbind, strsplit( sentrix, '_' ) ) )
  sentrix <- unique(cbind(i.link.list[,1],sentrix))
  sentrix <- sentrix[match(master.meth.id,sentrix[,1]),]
  sample.annotation <- data.frame(Sample_ID=as.character(sentrix[,1]),
                                  Sentrix_ID=as.character(sentrix[,2]),
                                  Sentrix_Position=as.character(sentrix[,3]) )
  write.csv(sample.annotation,paste(projects[i],"/",projects[i],"_sample_annotation.csv",sep=""),row.names=F)

  # File path
  sample.annotation <- file.path(paste(projects[i],"/",projects[i],"_sample_annotation.csv",sep=""))
  rnb.options(import.table.separator=",")

  # Report directory
  command <- paste("rm -fr ",projects[i],"/","RnBeads_normalization/",sep="")
  system(command)
  report.dir <- file.path(paste(projects[i],"/","RnBeads_normalization/",sep=""))

  # Vanilla parameters
  rnb.options(identifiers.column="Sample_ID")

  # Multiprocess
  num.cores <- 20
  parallel.setup(num.cores)
  
  #idat files
  idat.dir <- file.path("/root/TCGA/tcgaBiolink")
  data.source <- c(idat.dir, sample.annotation)
  result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)


#   485577 probes
#        of which: 482421 CpG, 3091 CpH, and 65 rs
#Region types:
#          138358 regions of type tiling
#           31033 regions of type genes
#           31195 regions of type promoters
#           26662 regions of type cpgislands
filt0 <- rnb.execute.na.removal(result$rnb.set, 0)$dataset
filt05 <- rnb.execute.na.removal(result$rnb.set, 0.05)$dataset
filt1 <- rnb.execute.na.removal(result$rnb.set, 1)$dataset
###
mf0 <- meth(filt0)
mf1 <- meth(filt1)
saveRDS(mf0,"mf0.rds")
saveRDS(mf1,"mf1.rds")

x = apply(mf0, 2, function(x) sum(is.na(x)) )


# norm
  rnb.set.norm <- rnb.execute.normalization(result$rnb.set, method="swan",bgcorr.method="methylumi.noob")
  save.rnb.set(rnb.set.norm,
               path=paste("/root/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/rnb.set.norm_withNormal.RData",sep=""))

  meth.norm<-meth(rnb.set.norm,row.names=T)

# meth batch correction
########################################################
########################################################
# extract control probe intensities(exluding negative control probes)
qc_data<-qc(rnb.set.example)
# Remove the 614 negative control probes from the intensity
control.annotation <- rnb.get.annotation("controls450")
non_neg_ctr_prob <- rownames(control.annotation[control.annotation$Target!="NEGATIVE",])
# add green intensity for non-negativeControl probes
qc_data_ctr_noNeg <- qc_data[[1]][rownames(qc_data[[1]]) %in% non_neg_ctr_prob,]
rownames(qc_data_ctr_noNeg) <- gsub("$","_green",rownames(qc_data_ctr_noNeg),perl=TRUE) #change name to avoid conflict
# add red intensity for non-negativeControl probes
qc_data_ctr_noNeg <- rbind(qc_data_ctr_noNeg,qc_data[[1]][rownames(qc_data[[2]]) %in% non_neg_ctr_prob,])
# remove rows with NAs
qc_data_ctr_noNeg <- qc_data_ctr_noNeg[complete.cases(qc_data_ctr_noNeg), ]

#PCA of control-probe intensities
pca <- prcomp(na.omit(t(qc_data_ctr_noNeg)))
ctrlprobes.scores = pca$x
colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
dim(ctrlprobes.scores)
phe=as.data.frame(ctrlprobes.scores)
          
lfla=as.formula('beta[i, ] ~  phe$PC1_cp + phe$PC2_cp + 
    phe$PC3_cp + phe$PC4_cp + phe$PC5_cp + phe$PC6_cp + phe$PC7_cp + 
    phe$PC8_cp + phe$PC9_cp + phe$PC10_cp + phe$PC11_cp + phe$PC12_cp + 
    phe$PC13_cp + phe$PC14_cp + phe$PC15_cp + phe$PC16_cp + phe$PC17_cp + 
    phe$PC18_cp + phe$PC19_cp + phe$PC20_cp + phe$PC21_cp + phe$PC22_cp + 
    phe$PC23_cp + phe$PC24_cp + phe$PC25_cp + phe$PC26_cp + phe$PC27_cp + 
    phe$PC28_cp + phe$PC29_cp + phe$PC30_cp')

# regression
samples=colnames(beta)
nvar = nrow(beta)
res=matrix(ncol=ncol(beta), nrow=nrow(beta))
rownames(res)=rownames(beta)
colnames(res)=colnames(beta)
for(i in 1:nvar) {
	tryCatch({fit = lm(lfla,na.action=na.exclude)}, error = function(error) {return(NA)})
	if(!exists("fit")){
		res[i,colnames(as.matrix(beta))] = rep(NA, length(colnames(as.matrix(beta))))
	}else{
		res[i,rownames(as.matrix(fit$residuals))] = fit$residuals
       		rm(fit)
	}		
}
          
# PCA on the resulting regression residuals (excluding markers with missing data) and include PC 1 to 5 as linear predictors in the final regression model.
pca <- prcomp(t(na.omit(res)))
pca.scores = pca$x[,1:30]
          
########################################################
########################################################
          
          
          
# Download package tarball from CRAN archive

url <- "https://cran.r-project.org/src/contrib/Archive/bigpca/bigpca_1.0.3.tar.gz"
pkgFile <- "bigpca_1.0.3.tar.gz"
download.file(url = url, destfile = pkgFile)

ERROR: dependencies ‘reader’, ‘NCmisc’, ‘bigmemory’, ‘biganalytics’, ‘bigmemory.sri’, ‘irlba’ are not available for package ‘bigpca’
install.packages(c("reader", "NCmisc", "bigmemory"))
install.packages(c("biganalytics", "bigmemory.sri", "irlba"))

          
# Install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
unlink(pkgFile)

library(bigpca)

beta<-meth(rnb.set.norm)
bmat2 <- as.big.matrix(beta)
## calculate PCA ##
result2 <- big.PCA(bmat2,thin=FALSE)
corrected <- PC.correct(result2,bmat2)
corrected2 <- PC.correct(result2,bmat2,n.cores=20)

          
          
# refFreeEWASP test
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
rnb.set.norm <- load.rnb.set("/root/TCGA/Rnbeads/KIRP/RnBeads_normalization/rnb.set.norm_withNormal.RData.zip")

annotFile <- paste0("/root/TCGA/Rnbeads/KIRP/KIRP_TP53_mutation_info_withNormal.csv") 
TUMOR = read.csv(annotFile,header=TRUE)
TUMOR = as.character(TUMOR$Variant_Classification)
TUMOR[TUMOR!="NORMAL"] = "TUMOR"
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Tumor = TUMOR)

clinical_path <- paste0("/root/TCGA/Rnbeads/KIRP/KIRP_clinical_info.csv") 
clinical_harmonized <- read.csv(clinical_path)
#covariates annotation
gender <- c( as.character(clinical_harmonized$GENDER), rep("NORMAL",sum(TUMOR=="NORMAL")) ) 
race <- c( as.character(clinical_harmonized$RACE), rep("NORMAL",sum(TUMOR=="NORMAL")) ) 
age <- c( as.numeric(clinical_harmonized$AGE), rep(NA,sum(TUMOR=="NORMAL")) ) 
rnb.set.norm@pheno = data.frame(rnb.set.norm@pheno, Gender = gender, Race = race, Age = age)


num.cores <- 20
parallel.setup(num.cores)

rnb.options(differential.site.test.method="refFreeEWAS")
dmc <- rnb.execute.computeDiffMeth(rnb.set.norm,pheno.cols=c("Tumor"))
          
 comparison <- get.comparisons(dmc)[1]
  dmc_table <-get.table(dmc, comparison, "sites", return.data.frame=TRUE)
  dmp_table <-get.table(dmc, comparison, "promoters", return.data.frame=TRUE)

          
          

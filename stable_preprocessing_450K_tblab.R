
# exodus
# idats
# scp -r -P 60057 root@172.18.149.78:/root/TCGA/tcgaBiolink/ /home/rtm/TCGA
# harmonized data from cbioportal
# scp -r -P 60057 root@172.18.149.78:/root/TCGA/panCancer_2018 /home/rtm/TCGA
# mkdir /home/rtm/TCGA/ex3
# mkdir /home/rtm/TCGA/ex3/Rnbeads

setwd("/home/rtm/TCGA/ex3/Rnbeads")
suppressMessages(library(RnBeads))
options(bitmapType="cairo")
options(scipen=999)
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
rnb.options(differential.site.test.method="refFreeEWAS")
library(TCGAbiolinks)

## preprocessing
#idat files
idat.dir <- file.path("/home/rtm/TCGA/tcgaBiolink/")

#### Link list processing ####
# master link list
link.list <- read.table("/home/rtm/TCGA/tcgaBiolink/idat_filename_case.txt",header=T,sep="\t")
# remove double associated methylation profiles
link.list <- link.list[grep(",",link.list$cases,invert=TRUE),]
# get project names
projects <- c("BLCA","BRCA","CHOL","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PRAD","THCA","UCEC")
# Make sure it matches the panCan data
panCan.dir <- list.dirs(path = "/home/rtm/TCGA/panCancer_2018", full.names = TRUE, recursive = FALSE)
for(i in 1:length(projects)){  if(grep(tolower(projects[i]),panCan.dir) > 0){print(paste(projects[i],": OK"))}  }

# Get all genes that can have mutations across all projects
master.gene.list <- read.table(pipe("cat /home/rtm/TCGA/panCancer_2018/*/data_mutations_mskcc.txt|tail -n +2|cut -f1|sort|uniq"),sep="\n")
master.gene.list <- as.character(master.gene.list[,1])
##############################
#### Project-based processing ####
for(i in 1:3){
  #project specific link.list
  print(paste("Parsing",projects[i], i,"Files"))
  
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
  print(paste("Starting 450K processing of",projects[i],i))

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
  num.cores <- 6
  parallel.setup(num.cores)
  
  #idat files
  idat.dir <- file.path("/home/rtm/TCGA/tcgaBiolink")
  data.source <- c(idat.dir, sample.annotation)
  result <- rnb.run.import(data.source=data.source,data.type="infinium.idat.dir", dir.reports=report.dir)

  # Annotate Pheno table 
  TUMOR = as.character(mut.file.p53$Variant_Classification)
  TUMOR[TUMOR!="NORMAL"] = "TUMOR"

  #covariates annotation
  gender <- as.character(c( as.character(clinical$GENDER), rep(NA,sum(TUMOR=="NORMAL")) ) )
  race <- as.character(c( as.character(clinical$RACE), rep(NA,sum(TUMOR=="NORMAL")) ) )
  age <- as.numeric(c( as.numeric(clinical$AGE), rep(NA,sum(TUMOR=="NORMAL")) ) ) 
	
  # normal sample annotation
  clinical <- GDCquery_clinic(project = paste0("TCGA-",projects[i]), type = "clinical")
  
  norm_samp<-rownames(rnb.set.norm@pheno)[TUMOR=="NORMAL"]

  parsed_clinical = data.frame(
      ID = clinical$bcr_patient_barcode,
      Gender = clinical$gender,
      Race = clinical$race,
      Age = clinical$age_at_index )
	
  meth.norm.id <- data.frame( do.call( rbind, strsplit( (norm_samp), '-' ) ) )
  meth.norm.id <- paste(meth.norm.id[,1],meth.norm.id[,2],meth.norm.id[,3],sep="-")

  xx_normal_clinical_data <- parsed_clinical[ match( meth.norm.id, as.character(parsed_clinical$ID) ) ,]

  age[which(TUMOR=="NORMAL")] = xx_normal_clinical_data$Age
  gender[which(TUMOR=="NORMAL")] = as.character(xx_normal_clinical_data$Gender)
  race[which(TUMOR=="NORMAL")] = as.character(xx_normal_clinical_data$Race)

  gender = tolower(gender)
  race = tolower(race)
	
  # adding annotation into rnb.set
  result$rnb.set@pheno = data.frame(result$rnb.set@pheno, Tumor = TUMOR, Gender = gender, Race = race, Age = age)

  print(paste("Finished with dataHarmonization of",projects[i], i))
  print(paste("Filtering probes callRates of",projects[i], i))
  # filter callRate on rows (methylation Probes) keeping rows with over 98% call rate
  rnb.set.rowCallRate <- rnb.execute.na.removal(result$rnb.set, 0.02)$dataset
	
  # free RAM 1
  rm(result)
  rm(mut.file)
  rm(mut.file.ix)
  rm(mut.mat)
  rm(mut.mat.ix)
  rm(CNA)
  rm(parsed_clinical)
  rm(clinical)

  # filter probe detection pvalue <= 10e-16
  pvals <- dpval(rnb.set.rowCallRate, row.names=TRUE)
  pvals <- pvals<=10e-16

  # filter samples with callrate <= 98%
  rnb.set.sampleRMV = remove.samples( rnb.set.rowCallRate, samples(rnb.set.rowCallRate)[ (colSums(pvals)/dim(pvals)[1])<.98 ] )

  # free RAM 2
  rm(rnb.set.rowCallRate)
  rm(pvals)

  print(paste("Starting normalization of:",projects[i], i))
  #Normalization and background correction
  rnb.set.norm <- rnb.execute.normalization(rnb.set.sampleRMV, method="swan",bgcorr.method="methylumi.noob")
	
  # free ram 3
  rm(rnb.set.sampleRMV)

  print(paste("Writting normalized file:",projects[i], i))
  # write rnb.set.norm
  save.rnb.set(rnb.set.norm,
               path=paste("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/rnb.set.norm_withNormal.RData",sep=""))

  # write beta
  meth.norm<-meth(rnb.set.norm,row.names=T)
  saveRDS(meth.norm, paste("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/betaVALUES_withNormal.rds",sep=""))
  print(paste("Samples and probes for ", i , projects[i], dim(meth.norm) ) )
	
  # write mval
  mval.norm <- mval(rnb.set.norm,row.names=T)
  saveRDS(mval.norm, paste("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/mVALUES_withNormal.rds",sep=""))
	
  # free ram 4
  rm(mval.norm)

  #### CPACOR
  print(paste("Starting CPACOR of:",projects[i],i))
  # extract control probe intensities(exluding negative control probes)
  qc_data<-qc(rnb.set.norm)
  # Remove the 614 negative control probes from the intensity
  control.annotation <- rnb.get.annotation("controls450")
  non_neg_ctr_prob <- rownames(control.annotation[control.annotation$Target!="NEGATIVE",])
  # add green intensity for non-negativeControl probes
  qc_data_ctr_noNeg <- qc_data[[1]][rownames(qc_data[[1]]) %in% non_neg_ctr_prob,]
  rownames(qc_data_ctr_noNeg) <- gsub("$","_green",rownames(qc_data_ctr_noNeg),perl=TRUE) #change name to avoid conflict
  # add red intensity for non-negativeControl probes
  qc_data_ctr_noNeg <- rbind(qc_data_ctr_noNeg, qc_data[[2]][rownames(qc_data[[2]]) %in% non_neg_ctr_prob,])
  # remove rows with NAs
  qc_data_ctr_noNeg <- qc_data_ctr_noNeg[complete.cases(qc_data_ctr_noNeg), ]

  #PCA of control-probe intensities
  pca <- prcomp(na.omit(t(qc_data_ctr_noNeg)))
  ctrlprobes.scores = pca$x
  colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
  phe=as.data.frame(ctrlprobes.scores)
  saveRDS(phe, paste("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/PC30_controProbeIntensity.rds",sep=""))
         
  lfla=as.formula('beta[i_pac, ] ~  phe$PC1_cp + phe$PC2_cp + 
    phe$PC3_cp + phe$PC4_cp + phe$PC5_cp + phe$PC6_cp + phe$PC7_cp + 
    phe$PC8_cp + phe$PC9_cp + phe$PC10_cp + phe$PC11_cp + phe$PC12_cp + 
    phe$PC13_cp + phe$PC14_cp + phe$PC15_cp + phe$PC16_cp + phe$PC17_cp + 
    phe$PC18_cp + phe$PC19_cp + phe$PC20_cp + phe$PC21_cp + phe$PC22_cp + 
    phe$PC23_cp + phe$PC24_cp + phe$PC25_cp + phe$PC26_cp + phe$PC27_cp + 
    phe$PC28_cp + phe$PC29_cp + phe$PC30_cp + rnb.set.norm@pheno$Age + rnb.set.norm@pheno$Gender + rnb.set.norm@pheno$Race ')

  print(paste("Doing CPACOR regresion of residual (with 30PCs:",projects[i],i))
  # regression
  beta = meth.norm
  samples=colnames(beta)
  nvar = nrow(beta)
  res=matrix(ncol=ncol(beta), nrow=nrow(beta))
  rownames(res)=rownames(beta)
  colnames(res)=colnames(beta)
  for(i_pac in 1:nvar) {
	tryCatch({fit = lm(lfla,na.action=na.exclude)}, error = function(error) {return(NA)})
	if(!exists("fit")){
		res[i_pac,colnames(as.matrix(beta))] = rep(NA, length(colnames(as.matrix(beta))))
	}else{
		res[i_pac,rownames(as.matrix(fit$residuals))] = fit$residuals
       		rm(fit)
	}		
  }
          
  print(paste("Doing PCA of the regressed residual (with 5PC):",projects[i],i))
  # PCA on the resulting regression residuals (excluding markers with missing data) and include PC 1 to 5 as linear predictors in the final regression model.
  pca <- prcomp(t(na.omit(res)))
  dim(pca$x)
  pca.scores = pca$x[,1:30]
  saveRDS(pca.scores, paste("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/","RnBeads_normalization/PC5_resultingResiduals.rds",sep=""))

  # add PCs
  rnb.set.norm@pheno <- data.frame(rnb.set.norm@pheno, phe[,1:30], pca.scores[,1:5])
	
  # Select Covariates
  rnb.options("covariate.adjustment.columns"=c("PC1_cp", "PC2_cp", "PC3_cp", "PC4_cp", "PC5_cp", "PC6_cp", "PC7_cp",
                                             "PC8_cp", "PC9_cp", "PC10_cp", "PC11_cp", "PC12_cp", "PC13_cp", "PC14_cp",
                                             "PC15_cp", "PC16_cp", "PC17_cp", "PC18_cp", "PC19_cp", "PC20_cp", "PC21_cp",
                                             "PC22_cp", "PC23_cp", "PC24_cp", "PC25_cp", "PC26_cp", "PC27_cp", "PC28_cp",
                                             "PC29_cp", "PC30_cp", "PC1", "PC2", "PC3", "PC4", "PC5", "Gender", "Race", "Age" ))

  # Run differential meth analysis
  print(paste("Running differential methAnlysis with refFreeEWAS:",projects[i],i))
  dmc <- rnb.execute.computeDiffMeth(rnb.set.norm,pheno.cols=c("Tumor"))
  # extract the vals        
  comparison <- get.comparisons(dmc)[1]
  dmc_table <-get.table(dmc, comparison, "sites", return.data.frame=TRUE)
  dmp_table <-get.table(dmc, comparison, "promoters", return.data.frame=TRUE)

  write.table(comparison,paste0("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/RnBeads_normalization/",projects[i],"_comparison.txt"))
  write.csv(dmc_table,paste0("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/RnBeads_normalization/",projects[i],"_dmc_table.csv"))
  write.csv(dmp_table,paste0("/home/rtm/TCGA/ex3/Rnbeads/",projects[i],"/RnBeads_normalization/",projects[i],"_dmp_table.csv"))

}


############################

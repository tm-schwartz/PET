#Read in Libraries ---------------------------

rm(list=ls())

require(dplyr)
require(purrr)
require(tidyr)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Out .csv argument must be supplied (input file).n", call.=FALSE)
}
outfnBase <- args[1]



df<-readRDS("/scratch/Hudson/resource/mapDatamerge/MAP_fm_d20220301_m20220301.rds")
#read in qei
df1 <-read.csv("/data/h_vmac/hudson/asl/e1_redo_test/asl_epoch1_02_06_2023.csv",colClasses=c(map.id="character",session="character"))
print('warning setting quality score to 1 for messed up runs temporarily, remove this step for final analysis')
df1['quality_score'] <- if_else(is.na(df1$quality_score),1,df1$quality_score)
df <- merge(df,df1,by.x =c("map.id","epoch","asl.session.id"),by.y=c("map.id","epoch","session"),all=T)


#read in sCoV and threshold by value
dfcov <- read.csv("/data/h_vmac/hudson/asl/V3_cbfScovE1.csv",colClasses=c(session="character"))
df <- merge(df,dfcov,by.x =c("asl.session.id"),by.y=c("session"),all=T)

#cyst for 163 I believe 
df <- dplyr::filter(df,map.id!="163")
#remove the one person with dementia
df <- dplyr::filter(df,diagnosis.factor!="Dementia")
#limit to epoch 1
df <- dplyr::filter(df,epoch==1)


cols <- c("map.id", "epoch","asl.session.id","apoe4pos.factor","apoe4count","apoe2count","cdr.factor","diagnosis.factor","age",
          "raceethnicity.factor","fsrp.minus.age.points","sex.factor","raceethnicity.factor",
          "education","asl.rest.usable.hct.factor","asl.scan.date",
          "ma.grey.matter","ma.total.intracranial.vol","np.executive.composite","quality_score","grey_matter_sCOV")


df <- dplyr::select(df,cols)

#skipping this for now because still likely usable
#next remove unusable
#df$asl.rest.usable.hct.binary <- if_else(df$asl.rest.usable.hct.factor=="Yes",1,0)

df_base <- df

#usableCols <- c("asl.rest.usable.hct.binary")

#loop through usable variables and remove if any not usable
#for (usableCol in usableCols){
#  print(usableCol)
#  df_base <- df_base[df_base[usableCol]==1,]
#}

#write all Epoch 1 usable data to .csv file
#write.csv(x = df_base,file = outfnAll,quote=F,row.names = F)


print('thresholding by quality score 0.7')
df_base <- dplyr::filter(df_base,quality_score>0.7)

print('thresholding by sCoV 0.5')
df_base <- dplyr::filter(df_base,grey_matter_sCOV<0.5)

#filter for data we want

cols <- c("map.id", "epoch","asl.session.id","apoe4pos.factor","diagnosis.factor","age",
          "raceethnicity.factor","fsrp.minus.age.points","sex.factor","raceethnicity.factor",
          "education","diagnosis.factor","np.executive.composite","quality_score","grey_matter_sCOV")

#"asl_rest_prg_precentral_gyrus_hr" removing prg covar for now
#roiNonPVC <- read.csv("/scratch/Hudson/ASL/09_04_2022/asl_epoch1_09_04_2022.csv",colClasses=c(map.id="character",session="character"))


#df_base <- merge(df_base,roiNonPVC,by.x =c("map.id","epoch","asl.session.id"),by.y=c("map.id","epoch","session"),all=T)

df_base <- dplyr::select(df_base,cols)





#usableCols <- c("pwv.usable")
#loop through usable variables and remove if any not usable
#for (usableCol in usableCols){
#  print(usableCol)
#  df_base <- df_base[df_base[usableCol]==1,]
#}

df_base <- tidyr::drop_na(df_base)




#drop bad CBF run 
df_base <- dplyr::filter(df_base,map.id!="034")

#remove duplicates
df_base <- df_base[!duplicated(df_base[c("map.id")]),]


df_merged <- df_base

df_merged$intercept<-1

predictor <- "np.executive.composite"

variables_to_factorize <- c("sex.factor", "apoe4pos.factor", "raceethnicity.factor", "diagnosis.factor")
text <- c("\nVariables to be factorized: ", variables_to_factorize)

# Select the variables which need to be mean centered. Add "factorized" over here to mean center all the variables
# that were factorized. You can also type each column in by hand otherwise. If you have interaction variables, that need to be mean centered, add them here by name 
# or just add "interaction" to pull the variable. 
variables_to_meancenter <- c("age","fsrp.minus.age.points","np.executive.composite", "factorized")
text <- c("\nVariables to be mean centered: ", variables_to_meancenter)

default_columns <- c("map.id", "epoch", "asl.session.id")


# First step to create factor variables for sex, apoe4, race/ethnicity, diagnosis
  #if there is a 3 or 4 level factor other than diagnosis, you will have to code that manually
text <- c("\nFactorizing variables: ", variables_to_factorize)
#write(paste(text, collapse=" "), summary_filename, append=TRUE)

for (var in variables_to_factorize) {
  if (var != "diagnosis.factor"){
  new_column_name <- paste(var, "factorized", sep = "_")
  most_common_in_column <- names(which.max(table(df_merged[[var]]))) 
  #write(sprintf("\nVariable: %s, %s code = 0", var, most_common_in_column), summary_filename, append=TRUE)
  df_merged[,new_column_name] <- ifelse(df_merged[[var]] == most_common_in_column, 0, 1)} 
  else {
    df_merged$mci_factorized <- ifelse(df_merged$diagnosis.factor=="MCI", 1, 0)
    df_merged$aar_factorized <- ifelse(df_merged$diagnosis.factor!="MCI" & df_merged$diagnosis.factor!="Normal", 1, 0)
  }
}

each_df <- df_merged


  for (var in variables_to_meancenter){
    col_name <- names(each_df)[grepl(sprintf("%s$",var),names(each_df))]
    new_var_name <- paste(col_name, "centered", sep='_')
    each_df[, new_var_name] <- scale(each_df[col_name], center=TRUE, scale=FALSE) # Mean center the columns in each dataframe
    #assign(list_of_dataframes[i], each_df) # Assign the updated data to the original dataframe
  }
  text <- c("Variables being mean centered:", variables_to_meancenter)
  #write(paste(text, collapse = " "), summary_filename, append=TRUE)

full_filename <- paste(outfnBase, "full_design.csv", sep="/")
write.csv(each_df, file=full_filename, quote = FALSE)

designmatrix<-select(each_df, map.id, asl.session.id, intercept, names(each_df)[grepl("centered",names(each_df))])
  # If interaction term is present, pull that down at the end. Else pull down the mean centered predictor.  
  if (exists("interactions")){
    interaction_centered <- sprintf("%s_centered",interactions)
    designmatrix <- designmatrix%>%select(-interaction_centered,everything())
  } else {
    predictor_centered <- sprintf("%s_centered", predictor)
    designmatrix <- designmatrix%>%select(-predictor_centered,everything())
  }
  csv_filename <- paste(outfnBase, "finaldesign.csv", sep="/")
  write.csv(designmatrix, file=csv_filename, quote = FALSE)
  #write(sprintf("\nWriting design matrix to CSV: %s", csv_filename), summary_filename, append=TRUE)

  # Write out the design matrix dataframe from the intercept column onwards 
  fsl_design_df <- designmatrix[,c(grep("intercept", names(designmatrix)):length(names(designmatrix)))]
  fsl_design_matrix <- paste(outfnBase, "fsl_design.mat", sep="/")
  fsl_contrast_matrix <- paste(outfnBase, "fsl_contrast.con", sep="/")
  # Write out some stuff that fsl needs at the beginning 
  write(sprintf("/NumWaves\t%d", ncol(fsl_design_df)),fsl_design_matrix)
  write(sprintf("/NumPoints \t%d", nrow(fsl_design_df)),fsl_design_matrix, append=TRUE)
  text <- c("/PPheights",rep(1,ncol(fsl_design_df)))
  write(paste(text, collapse = " "), fsl_design_matrix, append=TRUE)
  write("/Matrix",fsl_design_matrix, append=TRUE)
  write.table(fsl_design_df, fsl_design_matrix, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write out the contrast matrix
  write(sprintf("/NumWaves\t%d", ncol(fsl_design_df)),fsl_contrast_matrix)
  write(sprintf("/NumContrasts \t%d", 2),fsl_contrast_matrix, append=TRUE)
  text <- c("/PPheights",rep(1,ncol(fsl_design_df)))
  write(paste(text, collapse = " "), fsl_contrast_matrix, append=TRUE)
  write("/Matrix",fsl_contrast_matrix, append=TRUE)
  text <- c(rep(0,ncol(fsl_design_df)-1),1)
  write(paste(text, collapse = " "), fsl_contrast_matrix, append=TRUE)
  text <- c(rep(0,ncol(fsl_design_df)-1),-1)
  write(paste(text, collapse = " "), fsl_contrast_matrix, append=TRUE)

#write.csv(x = df_base,file = outfnAnalysis,quote=F,row.names = F)




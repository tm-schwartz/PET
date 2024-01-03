#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


#Load package
library(parallel)
library(RNifti)
library(tools)
library(FRACTION)
#Uncomment nlme library if running mixed effects
#library(nlme)


#Model function
runModel <- function(y,vw_reg,covDF) {
  #If all of y is NA, we skip estimation
  if (all(is.na(y))) {
    return(NA)
  }
  covDF$y <- y
  covDF$vw_reg <- vw_reg
  m <- -9999

  #Example lm with voxelwise-regressor
  try(m <- lm(y ~ vw_reg + age+fsrp.minus.age.points+education+sex.factor+raceethnicity.factor+apoe4pos.factor+pwv+diagnosis.factor, data=covDF,na.action=na.omit),silent=TRUE)

  #try(m <- lm(y ~ age, data=covDF,na.action=na.omit),silent=TRUE)

  #Example lme with voxelwise-regressor
  #cntrl<-lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "optim")
  #try(m <- lme(y ~ echo.ef.calc_base*nil_scan_date_diff_yr*diagnosis.factor_base + education_base + sex.factor_base + raceethnicity.factor_base + fsrp.minus.age.points_base + apoe4pos.factor_base + age_base + vw_reg, data = covDF, na.action=na.omit, random=~1+nil_scan_date_diff_yr | map.id, control = cntrl),silent=TRUE)

  return(m)
}


#This script is used to conduct voxel wise imaging analysis on a 4D imaging chunk.
#Usage chunkProc.R inputImg chunkFolder chunkNumber nChunks maskImg covariateFile outPath
#inputImg: input image for analysis
#chunkFolder: folder with image chunk
#vwRegImg: folder with voxelwise-regressor image for analysis
#vwRegChunkFolder: folder with voxelwise-regressor chunks.
#chunkNumber: number of chunk you want to process from 0 to # of chunks - 1
#nChunks: the total number of pieces your input image waas broken into
#maskImg: image to mask the input image.
#covariateFile: file containing covariate information
#outPath - path to put output images

#uncomment and specify args below if running within Rstudio
#args = c("/projDir/templates/Merged_Smoothed.nii.gz",
#          "/projDir/processingInR/merged_chunks",
#          "/projDir/templates/Merged_vwReg.nii.gz",
#          "/projDir/processingInR/vwReg_chunks",
#          0,
#          52,
#          "/projDir/templates/GM_Average_Thr.nii.gz",
#          "/projDir/DataOrganization/longitudinal_dataframe.csv",
#          "/projDir/processingInR/model_chunks")

parseCL <- function(args) {
  #Parse command line arguements 
  #Returns a vector containing the parsed arguements [1] = inputImg [2] = chunkFolder [3] = vwRegImg
  #[4] = vwRegChunkFolder [5] = chunkNumber [6] nChunks [7] = maskImg [8] = covariateFile, [9] = outPath
  # test if there are seven arguments: if not, return an error
  if (length(args)!=9) {
    stop("Usage chunkProc.R inputImg chunkFolder chunkNumber nChunks maskImg covariateFile outPath")
  } else {
    inImg = args[1]
    if (!file.exists(inImg)) {
      stop("Input image file does not exist")
    }
    chunkFolder = args[2]
    if (!dir.exists(chunkFolder)) {
      stop("Folder containing image slices does not exist")
    }
    vwRegImg = args[3]
    if (!file.exists(inImg)) {
      stop("Voxelwise-Regressor image file does not exist")
    }
    vwRegChunkFolder = args[4]
    if (!dir.exists(vwRegChunkFolder)) {
      stop("Folder containing voxelwise-regressor slices does not exist")
    }
    chunkNumber = as.numeric(args[5])
    if (!is.wholenumber(chunkNumber) | !chunkNumber >= 0) {
      stop("Chunk number refers to the numeric value assigned to your image chunk. This should be an integer greater than or equal to 0.")
    }
    nChunks = as.numeric(args[6])
    if (!is.wholenumber(nChunks) | !nChunks > 0) {
      stop("nChunks refers to the number of pieces your input image was broken into. This should be an integer greater than 0.")
    }
    maskImg = args[7]
    if (!file.exists(maskImg)) {
      stop("Mask image file does not exist")
    }
    covFile = args[8]
    if (!file.exists(covFile)) {
      stop("Covariate file does not exist")
    }
    outPath = args[9]
    if (!dir.exists(outPath)) {
      stop("Output path does not exist")
    }
  }
  parsedArgs <- c(inImg,chunkFolder,vwRegImg,vwRegChunkFolder,chunkNumber,nChunks,maskImg,covFile,outPath)
  return(parsedArgs)
} 


#Main Script

{
  parsedArgs<-parseCL(args)
  inImg <- parsedArgs[1]
  #Get Dimensions of file, 1st element is number of dimensions 
  inImgHeader <- niftiHeader(inImg)
  baseName <- basename(inImg)
  baseName <- file_path_sans_ext(baseName,compression = TRUE)
  nDim = inImgHeader$dim[1]
  inImgDim <- inImgHeader$dim[2:(nDim+1)]
  if (length(inImgDim) != 4) {
    stop("Only 4D input files supported")
  }
  #Get number of slices and volumes in 4d file. Assume z direction to be used in determining number of slices
  nSlices <- inImgDim[3]
  nVolumes <- inImgDim[4]
  
  chunkFolder <-  parsedArgs[2]

  vwRegImg <- parsedArgs[3]
  vwRegImgHeader <- niftiHeader(vwRegImg)
  vwRegImgbaseName <- basename(vwRegImg)
  vwRegImgbaseName <- file_path_sans_ext(vwRegImgbaseName,compression = TRUE)
  nDim = vwRegImgHeader$dim[1]
  vwRegImgDim <- vwRegImgHeader$dim[2:(nDim+1)]
  if (length(vwRegImgDim) != 4) {
    stop("Only 4D voxelwise-regressor files supported")
  }
  #The script will only check dimensions from headers in arguements. It will not compare dimensions of chunks. Assumes you have split your voxelwise-regressor and outcome image into an equal number and size of chunks.
  if (!all(inImgDim == vwRegImgDim)) {
    stop("Input image dimensions and voxelwise-regressor image dimensions do not match")
  }
  
  vwRegChunkFolder <-  parsedArgs[4]
  chunkNumber <- as.numeric(parsedArgs[5])
  nChunks <- as.numeric(parsedArgs[6])

  #Gather mask dimension information
  maskImg <- parsedArgs[7]
  maskImgHeader <- niftiHeader(maskImg)
  nDim <- maskImgHeader$dim[1]
  maskImgDim <- maskImgHeader$dim[2:(nDim+1)]
  if (length(maskImgDim) != 3) {
    stop("Only 3D mask files supported")
  }
  if(!all(inImgDim[1:3] == maskImgDim)) {
    stop("Mask does not have the same volume dimensions as a volume in the input image")
  }
  maskImg <- readNifti(maskImg,internal=FALSE,volumes=NULL)
  
  covFile <- parsedArgs[8]
  covFile <- read.csv(covFile,sep=",")
  #covFile$cdrBin <- covFile$cdr.factor > 0
  #covFile$cdrBin <- factor(covFile$cdrBin)
  #levels(covFile$cdrBin) <- c("Normal","MCI")
  
  outPath <- parsedArgs[9]
  
  #Read in chunk 
  chunkNumberStr = sprintf("%04d",chunkNumber)
  chunkFileName <- sprintf("%s/chunk_%s_%s.RDS", chunkFolder,chunkNumberStr,baseName)
  inImgChunk <- readRDS(file = chunkFileName, refhook = NULL)
  chunkDim <- dim(inImgChunk)
  chunkThickness <- chunkDim[3]
  
  #Read in voxelwise-regressor chunk
  vwRegChunkFileName <- sprintf("%s/chunk_%s_%s.RDS",vwRegChunkFolder,chunkNumberStr,vwRegImgbaseName)
  vwRegImgChunk <- readRDS(file = vwRegChunkFileName, refhook = NULL)

  #Determine which slices the chunk represents in the 4d file
  if (chunkNumber == (nChunks-1)) {
    sliceEnd <- nSlices
    sliceBegin <- sliceEnd - chunkThickness + 1
  } else {
    sliceBegin <- (chunkThickness * chunkNumber) + 1
    sliceEnd <- chunkThickness * (1+chunkNumber)
  }
  #Create a chunk of the mask image by indexing the calculated slice numbers
  maskImgChunk <- maskImg[,,sliceBegin:sliceEnd]

    
  print("Finished importing image and mask information")
  {
    #Convert 0's in mask chunk to NA's
    start.time <- Sys.time()
    maskImgChunk[maskImgChunk==0] <- NA
    end.time <- Sys.time()
    print(sprintf("Time elapsed converting 0's in mask image chunk to NAs: %f",end.time-start.time))
  }
  
  {
    #Mask all volumes of the input chunk by the mask chunk. Only mask main chunk and not voxelwise-regressor. Model skips are only based on the input image, and do not rely on the voxelwise-regressor.
    start.time <- Sys.time()
    for (j in 1:nVolumes) {
      inImgChunk[,,,j] <- inImgChunk[,,,j] * maskImgChunk
    }
    end.time <- Sys.time()
    print(sprintf("Time elapsed masking all volumes of the image chunk: %f",end.time-start.time))
  }
  
  {
    #Gather each voxels value into a matrix across all volumes
    #Each row contains values from a single voxel. Each column represents a certain volume.
    start.time <- Sys.time()
    voxVectors <- apply(inImgChunk,4,function(x) x)
    end.time <- Sys.time()
    print(sprintf("Time elapsed creating matrix of voxel values across all volumes: %f",end.time-start.time))
  }

  {
    #Gather each voxelwise regressor value into a matrix across all volumes
    #Each row contains values from a single voxel. Each column represents a certain volume
    start.time <- Sys.time()
    vwRegVectors <- apply(vwRegImgChunk,4,function(x) x)
    end.time <- Sys.time()
    print(sprintf("Time elapsed creating matrix of voxelwise-regressor values across all volumes: %f",end.time-start.time))
  }

  gc()
  
  print("Estimating model")
  
  {
    #Estimate a model for each voxel in the image chunk using a specified number of cores
    ncores <- 2
    start.time <- Sys.time()
    cl <- makeCluster(ncores,type="FORK")
    #The output will be a list where each element is the model object for a given voxel. 
    #The order is the same as the order of the rows in the voxVectors matrix. Which is just the linear index
    #order of the voxels in the image chunk in dimensions 1,2,3.
    result <- parLapply(cl,seq_len(nrow(voxVectors)), function(x) {
      return(runModel(voxVectors[x,],vwRegVectors[x,],covDF=covFile))
    })
    stopCluster(cl)
    end.time <- Sys.time()
    print(sprintf("Time elapsed estimating models for each voxel: %f",end.time-start.time))
  }
  
  print("Saving model estimate")
  {
    #Save model object list to outpath
    outFN <- sprintf("%s/model_estimate_chunk_%s.RDS", outPath,chunkNumberStr)
    saveRDS(result,file=outFN,ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    print(sprintf("Time elapsed saving model object list: %f",end.time-start.time))
  }
  print("Completed chunk processing")
  
}




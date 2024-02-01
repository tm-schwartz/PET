#!/usr/bin/env Rscript
#This script takes in a set of estimated model lists and reconstructs them into beta, tStat, and standardized residual
#images

#Inputs:
#args[1] = original 4d file
#args[2] = path to folder containing model estimate chunks
#args[3] = Number of chunks
#args[4] = residual (error) degrees of freedom
#args[5] = output path
args = commandArgs(trailingOnly=TRUE)

#Reconstruct model estimate chunks and write out desired images

library(RNifti)
library(tools)
library(FRACTION)

#uncomment and specify args below if running in Rstudio
#args = c("/projDir/templates/Merge_Smoothed.nii.gz",
#         "/projDir/processingInR/model_chunks",
#         52,
#         000,
#         "/projDir/processingInR/reconstructed_images")

parseCL <- function(args) {
  #Parse command line arguements 
  #Returns a vector containing the parsed arguements [1] = inputImg [2] = chunkFolder [3] = # of chunks
  #[4] residual (error) dof [5] = output path
  
  # test if there are 5 arguments: if not, return an error
  if (length(args)!=5) {
    stop("Usage: reconstructImgFromModelChunks.R inputImg chunkFolder nChunks outPath")
  } else {
    inImg = args[1]
    if (!file.exists(inImg)) {
      stop("Input image file does not exist")
    }
    chunkFolder = args[2]
    if (!dir.exists(chunkFolder)) {
      stop("Chunk folder does not exist")
    }
    nChunks = as.numeric(args[3])
    if (!is.wholenumber(nChunks) | !nChunks > 0) {
      stop("Number of chunks must be an integer greater than 0")
    }
    dof = as.numeric(args[4])
    if (!is.wholenumber(dof) | !dof > 0) {
      stop("Residual degrees of freedom must be an integer greater than 0")
    }
    outPath = args[5]
    if (!dir.exists(outPath)) {
      stop("Output path does not exist")
    }
  }
  parsedArgs <- c(inImg,chunkFolder,nChunks,dof,outPath)
  return(parsedArgs)
} 

#function to calculate standardized residuals
calcStdRes <- function(residualVector,dof) {
  #Residual sum of squares
  RSS <-  sum(residualVector * residualVector)
  #Residual mean squared error corrected for degrees of freedom
  RMS <-  RSS/dof
  stdResiduals <- residualVector / sqrt(RMS)
  return(stdResiduals)
}

#Parse command line
parsedArgs <- parseCL(args)
inImg <- parsedArgs[1]
chunkFolder <- parsedArgs[2]
nChunks <- as.numeric(parsedArgs[3])
dof <- as.numeric(parsedArgs[4])
outPath <- parsedArgs[5]


#Get 4d input header info 
inHeader <- niftiHeader(inImg)
#Get dimensions
nDim <- inHeader$dim[1]
if (nDim != 4 ) {
  stop("Image is not 4D")
}
inDim <- inHeader$dim[2:5]
nVoxInSlice <- inDim[1] * inDim[2]
nSlices <-inDim[3]
nVolumes <- inDim[4]
nVoxels <- nSlices * nVoxInSlice

#Build chunk filename
baseName <- basename(inImg)
baseName <- file_path_sans_ext(baseName,compression = TRUE)
nChunkVec <- seq(0,nChunks-1)
chunkFns <- sprintf("%s/model_estimate_chunk_%04d.RDS",chunkFolder,nChunkVec)
stdResChunkFns <- sprintf("%s/stdRes_array_chunk_%04d.nii.gz",outPath,nChunkVec)
#initialize outputs
betaVector <- rep(NA,nVoxels)
tStatsVector <- rep(NA,nVoxels)
pStatsVector <- rep(NA,nVoxels)
errVoxelVector <- rep(NA,nVoxels)

#Set variable indicating if summary file is written
finishedSummary = 0
#Set analytical summary text filename
analSummaryFn <- sprintf("%s/analyticalModelSummary.txt",outPath)

#curr slice will equal nSlices + 1 at end of loop
currSlice <- 1
totalErrors <- 0

print("Beginning image reconstruction from model chunk")
start.time <- Sys.time()
for (i in 1:nChunks) {
  #Determine where the chunk slices are in a given volume of the 4d input
  modelChunk <- readRDS(file = chunkFns[i],refhook = NULL)
  nVoxInChunk <- length(modelChunk)
  chunkThickness <- nVoxInChunk / nVoxInSlice
  endSlice <- currSlice + chunkThickness - 1
  beginIdx <- (currSlice * nVoxInSlice + 1) - nVoxInSlice
  endIdx <- beginIdx + nVoxInChunk - 1 
  print(sprintf("%s: Chunk: %d, Begin Slice: %d, End Slice: %d, Begin Voxel Linear Index: %d, End Linear Index: %d",Sys.time(),
                i,currSlice,endSlice,beginIdx,endIdx))
  
  #If summary file has not been created check if information can be gathered from a voxel within the given chunk
  if (finishedSummary == 0) {
    modelChunkLength <- sapply(modelChunk,length)
    #If a list element has a length > 1 it would indicate a model object
    if (any(modelChunkLength > 1)) {
      #Grab first xoxel meeting this criteria to pull info from
      exampleModel = modelChunk[[which(modelChunkLength > 1)[1]]]
      runFormula <- Reduce(paste,deparse(exampleModel$call))
      #Calculate a T-statistic height threshold based on a Z of 2.3 and specified degrees of freedom
      anal_height_threshold = qt(pnorm(2.3),df=dof)
      #Write information to a summary file
      write(runFormula,analSummaryFn)
      write(sprintf("Degrees of Freedom: %d", dof),analSummaryFn,append=TRUE)
      write(sprintf("Height Threshold: %6f", anal_height_threshold),analSummaryFn,append=TRUE)
      #Summary complete
      finishedSummary = 1
    }
  }  

  #Determine if there were any voxels that could not estimate a model and label them by a value of 1
  errVoxVectorChunk = sapply(modelChunk, function(x) { 
    len = length(x)
    #Count errors to be elements in model estimate list with a value of -9999 which aws assigned during chunk processing
    if (length(x) == 1 & all(!is.na(x))) {
      if (x == -9999) {
        return (1)
      }
    }
    return(0)
  })
   errVoxelVector[beginIdx:endIdx] <- errVoxVectorChunk
   print("Model Error Count:")
   print(sum(errVoxVectorChunk))
   totalErrors <- totalErrors + sum(errVoxVectorChunk)
   print("Total Error Count:")
   print(totalErrors)
  
  
  #List containing residuals for each voxel in chunk
  stdResListChunk <- lapply(modelChunk, function(x) {
    if ( all(is.na(x)) | (length(x) == 1 & x[[1]][1] == -9999) ) {
      return(rep(NA,nVolumes))
    } else {
      residuals <- residuals(x)

      ##WHR editing this part because it can get messed up if running a model where input 4d image contains NAs 
      #and are excluded from model (such as in LME)
      #we need to adjust for fact that resids returned can vary in length
      if (length(residuals) <nVolumes){
        nasInModel <- nVolumes - length(residuals)
        dof_adjusted <- dof - nasInModel
        positions_to_adjust <- which(is.na(x$data), arr.ind=TRUE)[,1]
        adjust <-T
      } else{
      dof_adjusted <- dof
      adjust <-F
      }

      stdResiduals <-calcStdRes(residuals,dof_adjusted)

      if (adjust){
        for (idx in positions_to_adjust){
          stdResiduals <- append(stdResiduals, NA, after=idx-1)
        }
      }

      return(stdResiduals)
    }
  })
  #Convert residual list to array with shape of the image chunk.
  #Unlisting needs to occur with array shape to keep all of the residuals of a voxel together
  stdResChunkArray <- array(unlist(stdResListChunk),c(nVolumes,inDim[1],inDim[2],chunkThickness))
  #Shift the first dimension to be the last to match image chunk dimensions.
  stdResChunkArray <- aperm(stdResChunkArray,c(2,3,4,1))
  writeNifti(stdResChunkArray, stdResChunkFns[i], template = inHeader, datatype = "double") 
  
  #Extract betas from chunk
  betaVectorOfChunk <- sapply(modelChunk, function(x) {
    if ( all(is.na(x)) | (length(x) == 1 & x[[1]][1] == -9999) ) {
      return(NA)
    } else {
      #Example extracting beta from lme model
      return(summary(x)$tTable["pre_post","Value"])

      #Example extracting beta from lm model
      #return(summary(x)$coefficients["pre_post:CognitivePositive","Estimate"])
    }
  })
  #assign to linear indicies of full betaVector for the final 3d image
  betaVector[beginIdx:endIdx] <- betaVectorOfChunk
  
  
  #Extract tStats from chunk
  tStatsVectorOfChunk <- sapply(modelChunk, function(x) {
     if ( all(is.na(x)) | (length(x) == 1 & x[[1]][1] == -9999) ) {
       return(NA)
     } else {

      #Example extracting t-statistic from lme model
      return(summary(x)$tTable["pre_post","t-value"])

      #Example extracting t-statistic from lm model
      #return(summary(x)$coefficients["pre_post:CognitivePositive","t value"])
     }
   })
  #assign to linear indicies of full tStatsVector for the final 3d image
  tStatsVector[beginIdx:endIdx] <- tStatsVectorOfChunk
  
  #Extract p values from chunk
  pStatsVectorOfChunk <- sapply(modelChunk, function(x) {
    if ( all(is.na(x)) | (length(x) == 1 & x[[1]][1] == -9999) ) {
      return(NA)
    } else {

      #Example extracting p-value from lme model
      return(summary(x)$tTable["pre_post","p-value"])

      #Example extracting p-value from lm model
      #return(summary(x)$coefficients["pre_post:CognitivePositive","Pr(>|t|)"])
    }
  })
  #assign to linear indicies of full tStatsVector for the final 3d image
  pStatsVector[beginIdx:endIdx] <- pStatsVectorOfChunk
  
  currSlice <- currSlice + chunkThickness
  #Ensure memory heavy objects below cleared including on 1st iteration
  rm(modelChunk,stdResListChunk,stdResChunkArray)
  #Free memory of objects above and memory with no reference pointers from previous iteration
  gc()
}

end.time <- Sys.time()
print(sprintf("Time elapsed reconstructing image statistics %s:",end.time-start.time))


#Output images:
#betas, tStats, pStats, and a mask of where a model could not be estimated 
print("Preparing and writing output images")
start.time <- Sys.time()

betaArray <- array(NA,inDim[1:3])
betaArray[,,] <- betaVector
betafn <- paste(outPath,"/betaWithNas.nii.gz",sep="")
writeNifti(betaArray, betafn, template = inHeader, datatype = "double")
betaArray[is.na(betaArray)] <-0
betafn <- paste(outPath,"/betaWithZeros.nii.gz",sep="")
writeNifti(betaArray, betafn, template = inHeader, datatype = "double")  

tStatsArray <- array(NA,inDim[1:3])
tStatsArray[,,] <- tStatsVector
tStatsfn <- paste(outPath,"/tStatWithNas.nii.gz",sep="")
writeNifti(tStatsArray, tStatsfn, template = inHeader, datatype = "double")     
tStatsArray[is.na(tStatsArray)] <-0
tStatsfn <- paste(outPath,"/tStatWithZeros.nii.gz",sep="")
writeNifti(tStatsArray, tStatsfn, template = inHeader, datatype = "double")     

pStatsArray <- array(NA,inDim[1:3])
pStatsArray[,,] <- pStatsVector
pStatsfn <- paste(outPath,"/pStatWithNas.nii.gz",sep="")
writeNifti(pStatsArray, pStatsfn, template = inHeader, datatype = "double")  
pStatsArray[is.na(pStatsArray)] <-0
pStatsfn <- paste(outPath,"/pStatWithZeros.nii.gz",sep="")
writeNifti(pStatsArray, pStatsfn, template = inHeader, datatype = "double")  

errVoxelArray <- array(NA,inDim[1:3])
errVoxelArray[,,] <- errVoxelVector
errVoxelArray[is.na(errVoxelArray)] <-0
errMaskfn <- paste(outPath,"/modelErrMask.nii.gz",sep="")
writeNifti(errVoxelArray, errMaskfn, template = inHeader, datatype = "uint8")

end.time <- Sys.time()
print(sprintf("Finished writing output images: %s",end.time-start.time))


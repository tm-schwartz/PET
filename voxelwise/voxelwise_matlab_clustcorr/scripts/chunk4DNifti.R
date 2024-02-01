#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#Split 4D Nifti File into smaller 4D chunks
#Input Args:
  #args[1] = input 4D file
  #args[2] = chunkThickness
  #args[2] = outPath: should be an existing folder

library(RNifti)
library(tools)
library(FRACTION)
#Uncomment if running in rstudio otherwise pass args as command line arguements
#args = c("/projDir/templates/Merged_Smoothed.nii.gz",
#         3,
#         "/projDir/processingInR/merged_chunks")

parseCL <- function(args) {
  #Parse command line arguements 
  #Returns a vector containing the parsed arguements [1] = inputImg [2] = chunkThickness [3] = outputPath
  
  # test if there are 3 arguments: if not, return an error
  if (length(args)!=3) {
    stop("Usage: chunk4DNifti.R inputImg chunkThickness outPath")
  } else {
    inImg = args[1]
    if (!file.exists(inImg)) {
      stop("Input image file does not exist")
    }
    chunkThickness = as.numeric(args[2])
    if (!is.wholenumber(chunkThickness) | !chunkThickness > 0) {
      stop("Chunk thickness refers to the number of slices in the output image chunks. This should be an integer greater than 0.")
    }
    outPath = args[3]
    if (!dir.exists(outPath)) {
      stop("Output path does not exist")
    }
  }
  parsedArgs <- c(inImg,chunkThickness,outPath)
  return(parsedArgs)
} 


{
  #Parse command line
  parsedArgs <- parseCL(args)
  inImg <- parsedArgs[1]
  chunkThickness <- as.numeric(parsedArgs[2])
  outPath <- parsedArgs[3]
  #Get input header info 
  inHeader <- niftiHeader(inImg)
  #Get dimensions, slice, and volume information
  nDim <- inHeader$dim[1]
  if (nDim != 4 ) {
    stop("Image is not 4D")
  }
  inDim <- inHeader$dim[2:5]
  nSlices <-inDim[3]
  nVolumes <- inDim[4]
  
  #Read in 4D file in segments of nVolumesPerRead
  nVolumesPerRead <- 100
  #Variable to keep track of amount left to read in loop
  nVolumesRem <- nVolumes
  #Put each segment as an element of list inImageSegments
  inImageSegments <- list()
  
  print("Reading in 4D Nifti...")
  start.time <- Sys.time()
  #Keep reading in nVolumes per read until no more volumes remain
  i <- 0
  while (nVolumesRem > 0 ) {
    if (nVolumesRem < nVolumesPerRead) {
      inImageSegments[[(i+1)]] <- readNifti(inImg,internal=FALSE,volumes = (nVolumesPerRead*i+1):(nVolumesPerRead*i+nVolumesRem))
      nVolumesRem <- nVolumesRem - nVolumesRem
    } else {
      inImageSegments[[(i+1)]] <- readNifti(inImg,internal=FALSE,volumes = (nVolumesPerRead*i+1):(nVolumesPerRead*i+nVolumesPerRead))
      nVolumesRem <- nVolumesRem - nVolumesPerRead
    }
    attr(inImageSegments[[(i+1)]],'.nifti_image_ptr') <- NULL
    gc()
    i <- i + 1
  }
  end.time <- Sys.time()
  timeDif <- difftime(end.time,start.time)
  print(sprintf("Time elapsed reading in file: %f %s",timeDif, units(timeDif)))
  
  #A chunk will be defined as a subset of voxels gathered across all volumes
  #Lets define chunks to be in the z direction for simplicity
  #Each chunk will have dimensions inDim[1] x inDim[2] x chunkThickness x nVolumes (aka inDim[4])
  
  
  #calculate how many chunks make up the image
  nChunks <- ceiling(nSlices / chunkThickness)
  nChunkVec <- seq(0,nChunks-1)
  #Gather output filenames
  #We will label filenames from 0 to nChunks - 1, which is how fslsplit tends to number split images
  baseName <- basename(inImg)
  baseName <- file_path_sans_ext(baseName,compression = TRUE)
  chunkFns <- sprintf("%s/chunk_%04d_%s.RDS",outPath,nChunkVec,baseName)
  
  #Create chunks of the image.
  #This is done by appropriately indexing slices across the image segments from read in.
  #Save the gathered image slices as an rds file in outpath
  slicesRemaining <- nSlices
  print("Creating and saving image chunks")
  start.time <- Sys.time()
  for (i in 0:(nChunks-1)) {
    volsRemaining <- nVolumes
    #Exit condition when there are less slices remaining than the actual chunk thickness we are at the end of the image
    if (slicesRemaining < chunkThickness) {
      imgChunk <- array(NA,c(inDim[1],inDim[2],slicesRemaining,inDim[4]))
      for (j in 0:(length(inImageSegments)-1)) {
        if (volsRemaining < nVolumesPerRead) {
          if (volsRemaining == 1) {
            #If only one volume remains index as a 3d array
            imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+volsRemaining)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+chunkThickness)]
          } else {
            imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+volsRemaining)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+chunkThickness),1:volsRemaining]
          }
          volsRemaining = volsRemaining - volsRemaining
        } else {
          imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+nVolumesPerRead)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+slicesRemaining),1:nVolumesPerRead]
          volsRemaining = volsRemaining - nVolumesPerRead
        }
      }
      slicesRemaining = slicesRemaining - slicesRemaining
    } else {
      imgChunk <- array(NA,c(inDim[1],inDim[2],chunkThickness,inDim[4]))
      #Loop thorugh the image segments and assign the appropriate slices of each volume to the image chunk
      for (j in 0:(length(inImageSegments)-1)) {
        if (volsRemaining < nVolumesPerRead) {
          if (volsRemaining == 1) {
            #If only one volume remains index as a 3d array
            imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+volsRemaining)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+chunkThickness)]
          } else {
            imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+volsRemaining)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+chunkThickness),1:volsRemaining]
          }
          volsRemaining = volsRemaining - volsRemaining
        } else {
          imgChunk[,,,(j*nVolumesPerRead+1):(j*nVolumesPerRead+nVolumesPerRead)] <- inImageSegments[[(j+1)]][,,(i*chunkThickness+1):(i*chunkThickness+chunkThickness),1:nVolumesPerRead]
          volsRemaining = volsRemaining - nVolumesPerRead
        }
      }
      slicesRemaining = slicesRemaining - chunkThickness
    }
    saveRDS(imgChunk,file=chunkFns[i+1], ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    rm(imgChunk)
    gc()
  }
  end.time <- Sys.time()
  timeDif <- difftime(end.time,start.time)
  print(sprintf("Time elapsed creating and saving image chunks: %f %s",timeDif,units(timeDif)))
}

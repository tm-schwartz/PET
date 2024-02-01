#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#Restructure 4d image slice chunks to 4d image volume chunks
#Input Args:
#args[1] = input 4D file
#args[2] = chunkFolder: folder contianing image chunks
#args[3] = nChunks: number of total chunks
#args[4] = nVolumesPerOutput
#args[5] = outPath: an existing folder for outputs
library(RNifti)
library(tools)
library(FRACTION)

#Uncomment if running in rstudio otherwise pass args as command line arguements
#args = c("/projDir/processingInR/Merge_Smoothed.nii.gz",
#         "/projDir/processingInR/reconstructed_images",
#         52,
#         25,
#         "/projDir/processingInR/restructured_images")

parseCL <- function(args) {
  #Parse command line arguements 
  #Returns a vector containing the parsed arguements [1] = inputImg [2] = chunkFolder 
  #[3] nChunks [4] nVolumesPerOutput [5] outPath
  # test if there are 5 arguments: if not, return an error
  if (length(args)!=5) {
    stop("Usage restructureSliceChunkToVolumeChunk.R inputImg chunkFolder nChunks outPath")
  } else {
    inImg = args[1]
    if (!file.exists(inImg)) {
      stop("Input image file does not exist")
    }
    chunkFolder = args[2]
    if (!dir.exists(chunkFolder)) {
      stop("Folder containing image slices does not exist")
    }
    nChunks = as.numeric(args[3])
    if (!is.wholenumber(nChunks) | !nChunks > 0) {
      stop("Number of chunks must be an integer greater than 0")
    }
    nVolumesPerOutput = as.numeric(args[4])
    if (!is.wholenumber(nVolumesPerOutput) | !nVolumesPerOutput > 0) {
      stop("Number of output volumes per chunk must be an integer greater than 0")
    }
    outPath = args[5]
    if (!dir.exists(outPath)) {
      stop("Output path does not exist")
    }
  }
  parsedArgs <- c(inImg,chunkFolder,nChunks,nVolumesPerOutput,outPath)
  return(parsedArgs)
} 

{
  #Parse command line
  parsedArgs <- parseCL(args)
  inImg <- parsedArgs[1]
  chunkFolder <- parsedArgs[2]
  nInChunks <- as.numeric(parsedArgs[3])
  nVolumesPerOutput <- as.numeric(parsedArgs[4])
  outPath <- parsedArgs[5]
  
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
  
  #Build chunk filenames
  inChunkVec <- seq(0,nInChunks-1)
  #Gather output filenames
  #We will label filenames from 0 to nChunks - 1, which is how fslsplit tends to number split images
  inChunkFns <- sprintf("%s/stdRes_array_chunk_%04d.nii.gz",chunkFolder,inChunkVec)
  
  
  #Gather outptu filenames
  nOutChunks <- ceiling(nVolumes / nVolumesPerOutput)
  outChunkVec <- seq(0,nOutChunks-1)
  outChunkFns <- sprintf("%s/volume_chunk_stdRes_%04d.nii.gz",outPath,outChunkVec)
  
  
  currVolume <- 0
  #Create output chunks
  for (i in 0:(nOutChunks-1)) {
    beginVol <- i*nVolumesPerOutput + 1
    currSlice <- 0
    #Determine which volumes the output chunk will represent
    if (i == (nOutChunks-1)) {
      endVol <- nVolumes
      remVols <- endVol - beginVol + 1
      outChunk <- array(NA,c(inDim[1],inDim[2],inDim[3],remVols))
      print(sprintf("%s: Output Chunk: %d, Begin Volume: %d, End Volume: %d",Sys.time(),
                    i+1,beginVol,endVol))
    } else {
      endVol <- i*nVolumesPerOutput + nVolumesPerOutput
      outChunk <- array(NA,c(inDim[1],inDim[2],inDim[3],nVolumesPerOutput))
      print(sprintf("%s: Output Chunk: %d, Begin Volume: %d, End Volume: %d",Sys.time(),
                    i+1,beginVol,endVol))
    }
    #Loop through the input chunks made up of slices in the z direction
    #  and read in the respective volumes needed for the output chunk
    for (j in 0:(nInChunks-1)) {
      #Read nifti and keep c pointer since we will gc after saving each output chunk.
      #Since gc will not occur until after all input chunks are read
      #Max memory required to read input data will equal
      #will = 2 * mem to store image of dimensions inDim1 x inDim2 x inDim3 x nVolumesPerOutput
      inChunk <- readNifti(inChunkFns[j+1],internal=FALSE,volumes = beginVol:endVol)
      chunkThickness <- dim(inChunk)[3]
      beginSlice <- currSlice + 1
      endSlice <- currSlice + chunkThickness
      print(sprintf("%s: Finished Reading Input Chunk: %d, Begin Slice: %d, End Slice: %d",Sys.time(),
                    j+1,beginSlice,endSlice))
      #Once completely filled outChunk will require the memory of storing an array
      #with dimensions of inDim1 x inDim2 xinDim3 x nVolumesPerOutput
      outChunk[,,beginSlice:endSlice,] <- inChunk
      currSlice <- endSlice
    }
    #Save output as .nii.gz
    #Writing the chunk will take an additional amount of memory approximately equal to
    #storing an image with dimensions of inDim1 x inDim2 x inDim3 x nVolumesPerOutput
    writeNifti(outChunk, outChunkFns[i+1], template = inHeader, datatype = "double")
    currVolume <- endVol
    #free memory
    rm(inChunk,outChunk)
    gc()
  }
}


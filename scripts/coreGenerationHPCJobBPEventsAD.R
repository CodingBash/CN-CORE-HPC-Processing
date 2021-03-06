#
# This script generates CORES from input files in BP units (instead of SNP/bin units).
# This script is not interactive and is meant to be ran on a HPC job from a shell script.
#

args <- commandArgs(trailingOnly = TRUE)

#
# Set script arguments
#
event <- "A" 
outputCsv <- "coreTable.csv"
outputObj <- "newCOREobj.rds"
if (length(args) == 1){
	event <- args[1]
} else if (length(args) == 2){
	event <- args[1]
	outputCsv <- args[2]	
} else if (length(args) == 3){
	event <- args[1]
	outputCsv <- args[2]
	outputObj <- args[3]
}

setwd("~/code/hN_core_artifacts/scripts")
source("coreGenerationLibrary.R")
source("helperFunctions.R")

# OVERWRITING EVENTS TO USE BOTH A and D
# TODO: Include both events via cmd line argument
events <- c("A", "D")
sample_dir <- "resources/FACETS_Reference_hN31_7_28_18_2/"
reference <- "hN31"

#
# Get CORE input
#
cd_core()
samples <- load_samples(classes = c("T", "F", "M"), sampleList = "./resources/sampleList.csv")
chromosomeSizes <- readRDS("./resources/chromosomeSizes.rds")

inputCORESegments <- selectSegmentsWithEvents(events = events, samples = samples, chromosomeSizes = chromosomeSizes, dir = sample_dir, sample_subdir="/", reference=reference, extension = "cnv.facets.v0.5.2.txt", inSampleFolder = TRUE, rescaleInput = TRUE, ampCall = 0.2, delCall = -0.235)


inputCOREBoundaries <- generateInputCOREBoundaries(chromosomeSizes)
print("Prepared all inputs - now running CORE")

#
# Run CORE
#
outputCOREobj <- runCORE(inputCORESegments, inputCOREBoundaries, distrib="Grid", maxmark=500, nshuffle=500, seedme=123, njobs=4)
print("CORE run complete")

#
# Save CORE object
# WARNING: The coreTable in the CORE object may not be in chromosomal location units. See CORE table output
#
saveRDS(outputCOREobj, outputObj)
print(paste("Saved CORE obj to", outputObj, sep = ""))


#
# Save CORE table
#
CORETable <- retrieveCORETable(outputCOREobj, chromosomeSizes, rescaleOutput = TRUE)
write.csv(CORETable, outputCsv)
print(paste("Saved coreTable as", outputCsv, sep = ""))

print("script complete")

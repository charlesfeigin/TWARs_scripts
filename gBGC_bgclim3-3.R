# Script to implement gBGC vs Null and gBGC+Sel vs gBGC LRTs
# Command line arguments are: 
#	1) the sequence alignment file in fasta format
#	2) the foreground branch of interest in the model file (e.g. Tcyn or Clup). This
#	   should match the name given in the neutral model file
#	3) the neutral model file
# Run script in a bash for loop over each fasta alignment of accelerated regions,
# using the bash redirect ">" to output results into an file. E.g.:
# for $file in *.fasta; do Rscript gBGC_bgclim3-3.R $file Tcyn mymodel.mod > $file.LRT 2> $file.err

# Load rphast library
library("rphast")

# Load command line arguments (see above)
args <- commandArgs(TRUE)

# Parse command line arguments
alnFasta <- args[1]
branch <- args[2]
mod <- args[3]

# Neutral model file hard-coded
model <- read.tm(mod)

# Read in the fasta alignment
align <- read.msa(alnFasta)

# Set model parameters and calculate likelihoods
modelData <- bgc.nucleotide.tests(align, model, branch, bgc.limits=c(3,3))

# Perform likelihood ratio tests
nullLikelihood <- modelData["null", "likelihood"]
gbgcLikelihood <- modelData["bgc", "likelihood"]
gbgcSelLikelihood <- modelData["sel+bgc", "likelihood"]
lrs <- c(gbgcLikelihood - nullLikelihood, + gbgcSelLikelihood - gbgcLikelihood)

# print results
print(paste(c(alnFasta, lrs), collapse = "    "))

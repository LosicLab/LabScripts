library(sciClone)
samples=list.files("vafs/")
system("mkdir -p scicloneOut_noCNV")

for ( myfile in samples[1:length(samples)] ) {
	if (file.exists(paste("scicloneOut_noCNV", myfile, sep="/"))) {
	} else {
	v1<-read.table(paste("vafs", myfile, sep="/"), header=T)
	try(writeClusterTable(sciClone(vafs=v1, sampleNames=myfile, minimumDepth=30), paste("scicloneOut_noCNV", myfile, sep="/")))
	}
}



#import library
library(circlize)

#read your bed files:
my_bed <- as.data.frame(read.table("evo_CDS.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed2 <- as.data.frame(read.table("evo_nCDS.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed3 <- as.data.frame(read.table("evo_intergenic.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
my_bed4 <- as.data.frame(read.table("evo_upstream.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))

# Init plot
circos.initializeWithIdeogram(my_bed4)

# Match the bed differences 
#bed_list = list(my_bed, my_bed2, my_bed3, my_bed4)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))

# specify each bed's methylation-density
circos.genomicDensity(my_bed, col = c("red"), track.height = 0.1)
circos.genomicDensity(my_bed2, col = c("blue"), track.height = 0.1)
circos.genomicDensity(my_bed3, col = c("purple"), track.height = 0.1)
circos.genomicDensity(my_bed4, col = c("orange"), track.height = 0.1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gridExtra))

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-bd", "--barcode_distribution_count", type="character", help="Counts of duplicates per barcode group", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-out", "--out_folder", type="character", help="Output folder", nargs=1, required=TRUE)


# reading parameters
args <- parser$parse_args()
group1 <- args$barcode_distribution_count
OUTF <- args$out_folder
print (OUTF)



# Ploting
setwd(OUTF)

GROUP1 <- as.data.frame(fread(group1,sep = '\t'))
colnames(GROUP1) <- c("COUNT")

GROUP_PLOT1 <- ggplot(GROUP1,aes(x=COUNT))+
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  ## scale_y_continuous(labels = percent_format()) #version 3.0.9
  scale_y_continuous(limits = c(0,1), breaks= seq(0,1,0.1)) +
  #scale_x_continuous(breaks = seq(1, 10, 1)) +
  #coord_cartesian(xlim=c(0, 10)) +
  ggtitle("Number of duplicates per barcode group") +
  theme_bw()+
  ylab("Fraction of barcode groups") +
  xlab("# Duplicates")


ggsave(file='DUPLICATES_PER_BARCODE_GROUP.pdf',plot=GROUP_PLOT1,useDingbats=FALSE)


# Some statistics

STATISTICS <- data.frame(mean(GROUP1$COUNT),median(GROUP1$COUNT))
colnames(STATISTICS) <- c("FAMILY(Mean)","FAMILY(Median)")

write.table(x = STATISTICS, file="BARCODE_STATS.txt",row.names = F,col.names = T,sep='\t', quote = F)


# Summary values

COUNT <- as.data.frame(table(GROUP1$COUNT))
colnames(COUNT) <- c("BARCODE_GROUP","COUNT")
#COUNT$Original_duplicates <- as.numeric(COUNT$BARCODE_GROUP) * as.numeric(COUNT$COUNT)

COUNT$Freq <- COUNT$COUNT/sum(COUNT$COUNT) 

write.table(x = COUNT, file="DUPLICATES_PER_BARCODE_GROUP.txt",row.names = F,col.names = T,sep='\t', quote = F)
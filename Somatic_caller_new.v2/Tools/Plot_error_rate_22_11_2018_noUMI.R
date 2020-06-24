suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library("seqinr")) 
suppressPackageStartupMessages(library("plyr")) 

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-bc1", "--bc1_file", type="character", help="TSV file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-out", "--out_folder", type="character", help="out_folder", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()

bc1 <- args$bc1_file

OUTF <- args$out_folder

setwd(OUTF)


#### SOME FUNCTIONS
ER_NT <- function(x){
  NT <- c("A","C","T","G")
  TAB <- NULL
  DP <- sum(as.numeric(x$DP))
  for (i in NT){
    X <- x[x$REF == i,]
    Z <- NT[!(NT == i)]
    for (y in Z){
      COUNT <- tryCatch(as.numeric(sum(X[y])), error = function(e) NA)
      TAB2 <- data.frame(i,y,sum(X$DP),COUNT,DP,COUNT/DP,COUNT/sum(X$DP))
      colnames(TAB2) <- c("REF","ALT","DP_REF","ALT_COUNT","DP_TOTAL", "ER_TOTAL","ER_REF")
      TAB <- rbind(TAB,TAB2)
    }
  }
  return(TAB)
}

NT_CHANGE <- function(x){
  x$Change <- apply(x,1,function(x) ifelse(x["REF"] %in% c("A","C"), paste(x["REF"],x["ALT"],sep='>'),paste(comp(x['REF'],forceToLower = F),comp(x['ALT'],forceToLower = F),sep='>')) )
  return(x)
}

LOADING <- function(DATA){
  NUCLEOTIDES <- c('A','C','T','G')
  
  DATA <- as.data.frame(fread(DATA))

  DATA$ALT_COUNT <- rowSums(DATA[,NUCLEOTIDES])-DATA$REFf-DATA$REFr
  DATA$DP <- rowSums(DATA[,NUCLEOTIDES])
  DATA$AB <- DATA$ALT_COUNT/DATA$DP_HQ

  #DATA <- DATA[DATA$AB < 0.05 & DATA$DP > 20,]
  DATA <- DATA[!(is.na(DATA$ALT_COUNT)),]

  return(DATA)
}

ERROR_RATE <- function(DATA){
  #DATA$AB <- DATA$ALT_COUNT/DATA$DP
  DATA <- DATA[DATA$DP > 10,]
  DATA <- DATA[!(is.na(DATA$ALT_COUNT)),]
  MEDIAN_DP <- median(DATA$DP)
  SUM_DP <- sum(as.numeric(DATA$DP))
  SUM_ALT <- sum(as.numeric(DATA$ALT_COUNT))    
  ER <- SUM_ALT/SUM_DP
  OUT <- data.frame(ER, nrow(DATA), SUM_ALT, SUM_DP)
  colnames(OUT) <- c('ER', 'Sites', 'ALT_COUNT', 'DP')
  return(OUT)
}

COLLAPSED_ER <- function(DATA){
  DATA[is.na(DATA)] <- 0
  COLLAPSED <- aggregate(cbind(ALT_COUNT,DP_REF) ~ Change, data = DATA, sum)
  COLLAPSED$ER <- COLLAPSED$ALT_COUNT/COLLAPSED$DP_REF
  return(COLLAPSED)
}

give.n <- function(x){
  return(c(y = 0, label = length(x)))
}


## Loading data
# Barcode 1x1
BC1 <- LOADING(bc1)

COV_PLOT <- ggplot(BC1, aes(x = DP_HQ))+
  geom_histogram(color = "dodgerblue3", fill = "lightskyblue2", binwidth = 1) +
  scale_fill_grey() +
  theme_bw() +
  labs(x='Depth of coverage', y='Count') +
  ggtitle("Depth of coverage")


## Finding potential germline calls
EXCLUDE <- BC1[BC1$AB >= 0.1,]

## Error Rates
# Barcode 1
BC1 <- BC1[!paste(BC1$CHROM, BC1$POS, sep = ';') %in% paste(EXCLUDE$CHROM, EXCLUDE$POS, sep = ';') & BC1$DP > 20 & BC1$AB < 0.1,]
ER1 <- ERROR_RATE(BC1)
colnames(ER1) <- c("ER", "Sites", "ALT_COUNT", "DP") 
ER1$Change <- 'All'
ER1$Type <- 'General error rate'

ER1_NT <- ER_NT(BC1)
ER1_NT <- NT_CHANGE(ER1_NT)
ER1_NT <- COLLAPSED_ER(ER1_NT)
colnames(ER1_NT) <- c("Change", "ALT_COUNT", "DP", "ER" )
ER1_NT$Type <- 'Nucleotide-specific error rate'

ERROR_RATES <- rbind.fill(ER1, ER1_NT)

### Getting error rates
ERROR_RATES$FILTER <- (ERROR_RATES$DP > 1e5)*1

ER_PASS <- ERROR_RATES[ERROR_RATES$FILTER == 1,]

COLORS <- c("All"= 'grey20', "A>C" = 'coral', "A>G" = 'darkgoldenrod1', "A>T" = 'darkseagreen', 
           "C>A" = 'deepskyblue2', "C>G" = 'azure4', "C>T" = 'tan4')

All_changes <- c("All", "A>C", "A>G", "A>T", 
                 "C>A", "C>G", "C>T")
PLOT_ER <- ggplot(ER_PASS,aes(x=Change,y=ER, fill = Change))+
  geom_bar(stat="identity") +
  #facet_wrap( ~ Type, ncol=2,  scales = "free_x", as.table = T) +
  scale_fill_manual(values=COLORS)+
  scale_x_discrete(limits=All_changes) + 
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Error rates")+
  xlab("Nucleotide change")+
  ylab("Error rate")



### PLOTING INFORMATION ERROR RATES
ggsave(file="Coverage.pdf",plot=COV_PLOT, useDingbats=FALSE)
ggsave(file='Error_rates.pdf',plot=PLOT_ER, useDingbats=FALSE)
write.table(ERROR_RATES,"Error_rates.txt",row.names = F,col.names = T,sep='\t', quote = F)



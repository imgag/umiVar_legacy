setwd("W:/users/ahmuyaf1/scripts/Somatic_caller_new.v2/Test")

TD1 <- 'DEDUPLICATED_DP1.tsv'
TD2 <- 'DEDUPLICATED_DP2.tsv'
TD3 <- 'DEDUPLICATED_DP3.tsv'
TD4 <- 'DEDUPLICATED_DP4.tsv'
parameters <- 'Parameters.txt'






## Indel genotype
ALL <- rbind(CALLING1, CALLING2, CALLING3, CALLING4)
ALL$COUNT <- 1

#ALL_del <- ALL[ALL$DEL > 0,]
if (nrow(ALL_del) > 0){
  Del <- aggregate(cbind(DEL, COUNT) ~ CHROM + POS + REF + DELg, ALL_del, sum)
  
  # first decide for the one more represented in more duplicate groups
  Del_max_count <- aggregate(COUNT ~ CHROM + POS, Del, max)
  colnames(Del_max_count) <- c('CHROM','POS','Max_count')
  
  Del <- merge(Del, Del_max_count, by = c('CHROM', 'POS'))
  Del <- Del[Del$COUNT == Del$Max_count,]
  
  # second choose for more represented variant
  Del_max_del <- aggregate(DEL ~ CHROM + POS, Del, max)
  colnames(Del_max_del) <- c('CHROM','POS','Max_del')
  
  Del <- merge(Del, Del_max_del, by = c('CHROM', 'POS'))
  Del <- Del[Del$DEL == Del$Max_del,]  
  
  Del <- Del[,c('CHROM', 'POS','DELg')]
  colnames(Del) <- c('CHROM', 'POS','Genotype')
  Del$Type <- 'DEL'
  
} else {
  Del <- numeric(0)
}

ALL_ins <- ALL[ALL$INS > 0,]
if (nrow(ALL_ins) > 0){
  Ins <- aggregate(cbind(INS, COUNT) ~ CHROM + POS + REF + INSg, ALL_ins, sum)
  
  # first decide for the one more represented in more duplicate groups
  Ins_max_count <- aggregate(COUNT ~ CHROM + POS, Ins, max)
  colnames(Ins_max_count) <- c('CHROM','POS','Max_count')
  
  Ins <- merge(Ins, Ins_max_count, by = c('CHROM', 'POS'))
  Ins <- Ins[Ins$COUNT == Ins$Max_count,]
  
  # second choose for more represented variant
  Ins_max_Ins <- aggregate(INS ~ CHROM + POS, Ins, max)
  colnames(Ins_max_Ins) <- c('CHROM','POS','Max_Ins')
  
  Ins <- merge(Ins, Ins_max_Ins, by = c('CHROM', 'POS'))
  Ins <- Ins[Ins$INS == Ins$Max_Ins,]  
  
  Ins <- Ins[,c('CHROM', 'POS','INSg')]
  colnames(Ins) <- c('CHROM', 'POS','Genotype')
  Ins$Type <- 'INS'
  
} else {
  Ins <- numeric(0)
}

INDEL <- rbind(Ins, Del)

LINE <- VARIANTS[14,]



for (i in row.names(VARIANTS)){
  
  J <- CALL(VARIANTS[i,], INDEL)
  print (c(i,J))
}

LINE <- Z

CALL <- function(LINE, INDEL, num_sites){
  LAT <- c("A", "C", "T", "G", "INS", "DEL")
  P <- c("PA", "PC", "PT", "PG", "PINS", "PDEL", "PINSo", "PDELo")
  P_adj <- c("PA_adj", "PC_adj", "PT_adj", "PG_adj", "PINS_adj", "PDEL_adj")
  FOR <- c("Af", "Cf", "Tf", "Gf", "INSf", "DELf")
  REV <- c("Ar", "Cr", "Tr", "Gr", "INSr", "DELr")
  ALL_id <- seq(1,length(P))
  
  num_sites <- max(as.numeric(num_sites), 0)
  
  REF <- as.character(LINE['REF'])
  
  ALL_ALT_COUNT <- c('A','C','T','G', 'INS', 'DEL', 'INSo','DELo')
  
  #LINE <- VARIANTS[1,]
  
  P_ADJ <- as.numeric(LINE[P_adj])
  
  id <- which(P_ADJ <= 0.05)
  
  #print (id)
  
  if (length(id) > 0){
    call <- NULL
    for (i in id){
      ALT <- LAT[i]
      
      if (ALT %in% c('DEL', 'INS')){
        
        GENOTYPE <- INDEL[paste(LINE['CHROM'], as.numeric(LINE['POS']), sep = '-') == paste(INDEL$CHROM, INDEL$POS, sep = '-') & 
                            INDEL$Type == ALT,]$Genotype
      } else {
        GENOTYPE <- paste(REF, ALT, sep = '>')
      }
      
      Ps <- as.numeric(LINE[P])
      FORWARD <- as.numeric(LINE[FOR[i]])
      REVERSE <- as.numeric(LINE[REV[i]])
      ALT_COUNT <- as.numeric(FORWARD+REVERSE)
      ALT_COUNT_all <- as.numeric(LINE[ALL_ALT_COUNT])
      ALT_COUNT_padj <- round(as.numeric(P_ADJ[i]),5)
      ALT_COUNT_p <- round(as.numeric(Ps[i]),5)
      
      AB <- round(ALT_COUNT/(as.numeric(LINE['REFf']) + as.numeric(LINE['REFr'])), 6)
      
      
      id_o <- ALL_id[ALL_id != id]
      Ps_out <- Ps[ALL_id[ALL_id != id]]
      OUT <- (prod(Ps_out, na.rm = T))
      OUT_adj <- round(p.adjust(OUT, n = num_sites, method = 'fdr'), 5)
      ALT_COUNT_o <- as.numeric(sum(ALT_COUNT_all[id_o]) - as.numeric(LINE['REFf']) - as.numeric(LINE['REFr']))
      
      
      call_temp <- paste(GENOTYPE, ALT_COUNT, paste(FORWARD, REVERSE,sep = '-'), AB, ALT_COUNT_p, ALT_COUNT_padj, ALT_COUNT_o, OUT_adj, sep = ';')
      
      call <- c(call, call_temp)
    }
    call <- paste(call, collapse = '|')
    
  } else {
    call <- 'No_call'
  }
  return(call)
}


Z <- VARIANTS
Z$CALL <- apply(Z, 1, function(x) CALL(x, INDEL, num_sites = num_sites))
Z$DP2 <- rowSums(Z[,c('A','C','T','G', 'INS', 'DEL', 'INSo','DELo')])
Z$DIFF <- Z$DP_HQ-Z$DP2
View(Z[,c('CHROM', 'POS', 'DP_HQ', 'DP2', 'INS','INSo','DEL','DELo', 'DIFF')])

FINAL <- Z[,c('CHROM','POS','REF','DP','DP_HQ', 'REFf', 'REFr','CALL')]

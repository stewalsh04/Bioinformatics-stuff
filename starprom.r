---
title: "Starprom"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

##Function to get reverse complement
```{r}
rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}
```

##Read in Promoter Data and barcode requences
```{r}
mydata <- read.csv(file="/home/steve/Starprom/all_joined.csv",sep="\t",header=FALSE)
mydata2 <- as.matrix(mydata)
barcodes <- read.csv(file="/home/steve/Starprom/barcodes.csv",sep="\t", header=TRUE)
promoters <- read.csv(file="/home/steve/Starprom/motif_promoters.csv",sep=",", header=TRUE)
promoters$promoter <- paste("GCCA", promoters$promoter, sep="")
promoters$promoter <- paste(promoters$promoter, "GAAG", sep="")
barcodes <- as.data.frame(sort(barcodes$Barcode))
colnames(barcodes) <- ("Barcode")
barcodes <- unique(barcodes)
reversecompbarcodes <- barcodes
barcodes$Barcode <- rc(barcodes$Barcode)  ##Get reverse comp of barcodes
```

##Add Column names and cut out all columns with NA ie no reads
```{r}
colnames(mydata) <- c("2021","2021","2021","2022","2022","2022","2041","2041","2041","2042","2042","2042","2081","2081","2081","2081","2082","2082","6811","6811","6811","6812","6812","6812","DMSO11","DMSO11","DMSO11","DMSO12","DMSO12","DMSO12","DMSO21","DMSO21","DMSO21","DMSO22","DMSO22","DMSO22","DMSO41","DMSO41","DMSO41","DMSO42","DMSO42","DMSO42","DMSO81","DMSO81","DMSO81","DMSO82","DMSO82","DMSO82","DMSO01","DMSO01","DMSO01","DMSO02","DMSO02","DMSO02","6821","6821","6821","6822","6822","6822","6841","6841","6841","6842","6842","6842","6881","6881","6881","6882","6882","6882")
#mydata <- mydata[complete.cases(mydata), ] #Sifts out rows with NA for every sample
```

##Split each Sample, select only barcodes in Starprom list, and reverse comp barcodes in sample counts
```{r}
data_2021 <- mydata[,c(3,2)]
colnames(data_2021) <- c("Barcode","Count_2021")
data_2021 <- merge(x = barcodes, y = data_2021[ , c("Barcode", "Count_2021")], by = "Barcode", all.x = TRUE)
data_2021$Barcode <- rc(data_2021$Barcode)
data_2021 <- data_2021[order(data_2021$Barcode),]
data_2021$Count_2021 <- as.numeric(as.character(data_2021$Count_2021))
#data_2021 <- aggregate(data_2021$data_2021, by = list(data_2021$Barcode), FUN=sum)

data_2022 <- mydata[,c(6,5)]
colnames(data_2022) <- c("Barcode","Count_2022")
data_2022 <- merge(x = barcodes, y = data_2022[ , c("Barcode", "Count_2022")], by = "Barcode", all.x = TRUE)
data_2022$Barcode <- rc(data_2022$Barcode)
data_2022 <- data_2022[order(data_2022$Barcode),]
data_2022$Count_2022 <- as.numeric(as.character(data_2022$Count_2022))
#data_2022 <- aggregate(data_2022$data_2022, by = list(data_2022$Barcode), FUN=sum)

data_2041 <- mydata[,c(9,8)]
colnames(data_2041) <- c("Barcode","Count_2041")
data_2041 <- merge(x = barcodes, y = data_2041[ , c("Barcode", "Count_2041")], by = "Barcode", all.x = TRUE)
data_2041$Barcode <- rc(data_2041$Barcode)
data_2041 <- data_2041[order(data_2041$Barcode),]
data_2041$Count_2041 <- as.numeric(as.character(data_2041$Count_2041))
#data_2041 <- aggregate(data_2041$data_2041, by = list(data_2041$Barcode), FUN=sum)

data_2042 <- mydata[,c(12,11)]
colnames(data_2042) <- c("Barcode","Count_2042")
data_2042 <- merge(x = barcodes, y = data_2042[ , c("Barcode", "Count_2042")], by = "Barcode", all.x = TRUE)
data_2042$Barcode <- rc(data_2042$Barcode)
data_2042 <- data_2042[order(data_2042$Barcode),]
data_2042$Count_2042 <- as.numeric(as.character(data_2042$Count_2042))
#data_2042 <- aggregate(data_2042$data_2042, by = list(data_2042$Barcode), FUN=sum)

data_2081 <- mydata[,c(15,14)]
colnames(data_2081) <- c("Barcode","Count_2081")
data_2081 <- merge(x = barcodes, y = data_2081[ , c("Barcode", "Count_2081")], by = "Barcode", all.x = TRUE)
data_2081$Barcode <- rc(data_2081$Barcode)
data_2081 <- data_2081[order(data_2081$Barcode),]
data_2081$Count_2081 <- as.numeric(as.character(data_2081$Count_2081))
#data_2081 <- aggregate(data_2081$data_2081, by = list(data_2081$Barcode), FUN=sum)

data_2082 <- mydata[,c(18,17)]
colnames(data_2082) <- c("Barcode","Count_2082")
data_2082 <- merge(x = barcodes, y = data_2082[ , c("Barcode", "Count_2082")], by = "Barcode", all.x = TRUE)
data_2082$Barcode <- rc(data_2082$Barcode)
data_2082 <- data_2082[order(data_2082$Barcode),]
data_2082$Count_2082 <- as.numeric(as.character(data_2082$Count_2082))
data_2082 <- aggregate(data_2082$Count_2082, by = list(data_2082$Barcode), FUN=sum)
colnames(data_2082) <- c("Barcode","Count_2082")

data_6811 <- mydata[,c(21,20)]
colnames(data_6811) <- c("Barcode","Count_6811")
data_6811 <- merge(x = barcodes, y = data_6811[ , c("Barcode", "Count_6811")], by = "Barcode", all.x = TRUE)
data_6811$Barcode <- rc(data_6811$Barcode)
data_6811 <- data_6811[order(data_6811$Barcode),]
data_6811$Count_6811 <- as.numeric(as.character(data_6811$Count_6811))
#data_6811 <- aggregate(data_6811$data_6811, by = list(data_6811$Barcode), FUN=sum)

data_6812 <- mydata[,c(24,23)]
colnames(data_6812) <- c("Barcode","Count_6812")
data_6812 <- merge(x = barcodes, y = data_6812[ , c("Barcode", "Count_6812")], by = "Barcode", all.x = TRUE)
data_6812$Barcode <- rc(data_6812$Barcode)
data_6812 <- data_6812[order(data_6812$Barcode),]
data_6812$Count_6812 <- as.numeric(as.character(data_6812$Count_6812))
data_6812 <- aggregate(data_6812$Count_6812, by = list(data_6812$Barcode), FUN=sum)
colnames(data_6812) <- c("Barcode","Count_6812")

data_DMSO11 <- mydata[,c(27,26)]
colnames(data_DMSO11) <- c("Barcode","Count_DMSO11")
data_DMSO11 <- merge(x = barcodes, y = data_DMSO11[ , c("Barcode", "Count_DMSO11")], by = "Barcode", all.x = TRUE)
data_DMSO11$Barcode <- rc(data_DMSO11$Barcode)
data_DMSO11 <- data_DMSO11[order(data_DMSO11$Barcode),]
data_DMSO11$Count_DMSO11 <- as.numeric(as.character(data_DMSO11$Count_DMSO11))
data_DMSO11 <- aggregate(data_DMSO11$Count_DMSO11, by = list(data_DMSO11$Barcode), FUN=sum)
colnames(data_DMSO11) <- c("Barcode","Count_DMSO11")

data_DMSO12 <- mydata[,c(30,29)]
colnames(data_DMSO12) <- c("Barcode","Count_DMSO12")
data_DMSO12 <- merge(x = barcodes, y = data_DMSO12[ , c("Barcode", "Count_DMSO12")], by = "Barcode", all.x = TRUE)
data_DMSO12$Barcode <- rc(data_DMSO12$Barcode)
data_DMSO12 <- data_DMSO12[order(data_DMSO12$Barcode),]
data_DMSO12$Count_DMSO12 <- as.numeric(as.character(data_DMSO12$Count_DMSO12))
data_DMSO12 <- aggregate(data_DMSO12$Count_DMSO12, by = list(data_DMSO12$Barcode), FUN=sum)
colnames(data_DMSO12) <- c("Barcode","Count_DMSO12")

data_DMSO21 <- mydata[,c(33,32)]
colnames(data_DMSO21) <- c("Barcode","Count_DMSO21")
data_DMSO21 <- merge(x = barcodes, y = data_DMSO21[ , c("Barcode", "Count_DMSO21")], by = "Barcode", all.x = TRUE)
data_DMSO21$Barcode <- rc(data_DMSO21$Barcode)
data_DMSO21 <- data_DMSO21[order(data_DMSO21$Barcode),]
data_DMSO21$Count_DMSO21 <- as.numeric(as.character(data_DMSO21$Count_DMSO21))
data_DMSO21 <- aggregate(data_DMSO21$Count_DMSO21, by = list(data_DMSO21$Barcode), FUN=sum)
colnames(data_DMSO21) <- c("Barcode","Count_DMSO21")

data_DMSO22 <- mydata[,c(36,35)]
colnames(data_DMSO22) <- c("Barcode","Count_DMSO22")
data_DMSO22 <- merge(x = barcodes, y = data_DMSO22[ , c("Barcode", "Count_DMSO22")], by = "Barcode", all.x = TRUE)
data_DMSO22$Barcode <- rc(data_DMSO22$Barcode)
data_DMSO22 <- data_DMSO22[order(data_DMSO22$Barcode),]
data_DMSO22$Count_DMSO22 <- as.numeric(as.character(data_DMSO22$Count_DMSO22))
data_DMSO22 <- aggregate(data_DMSO22$Count_DMSO22, by = list(data_DMSO22$Barcode), FUN=sum)
colnames(data_DMSO22) <- c("Barcode","Count_DMSO22")

data_DMSO41 <- mydata[,c(39,38)]
colnames(data_DMSO41) <- c("Barcode","Count_DMSO41")
data_DMSO41 <- merge(x = barcodes, y = data_DMSO41[ , c("Barcode", "Count_DMSO41")], by = "Barcode", all.x = TRUE)
data_DMSO41$Barcode <- rc(data_DMSO41$Barcode)
data_DMSO41 <- data_DMSO41[order(data_DMSO41$Barcode),]
data_DMSO41$Count_DMSO41 <- as.numeric(as.character(data_DMSO41$Count_DMSO41))
data_DMSO41 <- aggregate(data_DMSO41$Count_DMSO41, by = list(data_DMSO41$Barcode), FUN=sum)
colnames(data_DMSO41) <- c("Barcode","Count_DMSO41")

data_DMSO42 <- mydata[,c(42,41)]
colnames(data_DMSO42) <- c("Barcode","Count_DMSO42")
data_DMSO42 <- merge(x = barcodes, y = data_DMSO42[ , c("Barcode", "Count_DMSO42")], by = "Barcode", all.x = TRUE)
data_DMSO42$Barcode <- rc(data_DMSO42$Barcode)
data_DMSO42 <- data_DMSO42[order(data_DMSO42$Barcode),]
data_DMSO42$Count_DMSO42 <- as.numeric(as.character(data_DMSO42$Count_DMSO42))
data_DMSO42 <- aggregate(data_DMSO42$Count_DMSO42, by = list(data_DMSO42$Barcode), FUN=sum)
colnames(data_DMSO42) <- c("Barcode","Count_DMSO42")

data_DMSO81 <- mydata[,c(45,44)]
colnames(data_DMSO81) <- c("Barcode","Count_DMSO81")
data_DMSO81 <- merge(x = barcodes, y = data_DMSO81[ , c("Barcode", "Count_DMSO81")], by = "Barcode", all.x = TRUE)
data_DMSO81$Barcode <- rc(data_DMSO81$Barcode)
data_DMSO81 <- data_DMSO81[order(data_DMSO81$Barcode),]
data_DMSO81$Count_DMSO81 <- as.numeric(as.character(data_DMSO81$Count_DMSO81))
data_DMSO81 <- aggregate(data_DMSO81$Count_DMSO81, by = list(data_DMSO81$Barcode), FUN=sum)
colnames(data_DMSO81) <- c("Barcode","Count_DMSO81")

data_DMSO82 <- mydata[,c(48,47)]
colnames(data_DMSO82) <- c("Barcode","Count_DMSO82")
data_DMSO82 <- merge(x = barcodes, y = data_DMSO82[ , c("Barcode", "Count_DMSO82")], by = "Barcode", all.x = TRUE)
data_DMSO82$Barcode <- rc(data_DMSO82$Barcode)
data_DMSO82 <- data_DMSO82[order(data_DMSO82$Barcode),]
data_DMSO82$Count_DMSO82 <- as.numeric(as.character(data_DMSO82$Count_DMSO82))
data_DMSO82 <- aggregate(data_DMSO82$Count_DMSO82, by = list(data_DMSO82$Barcode), FUN=sum)
colnames(data_DMSO82) <- c("Barcode","Count_DMSO82")

data_DMSO01 <- mydata[,c(51,50)]
colnames(data_DMSO01) <- c("Barcode","Count_DMSO01")
data_DMSO01 <- merge(x = barcodes, y = data_DMSO01[ , c("Barcode", "Count_DMSO01")], by = "Barcode", all.x = TRUE)
data_DMSO01$Barcode <- rc(data_DMSO01$Barcode)
data_DMSO01 <- data_DMSO01[order(data_DMSO01$Barcode),]
data_DMSO01$Count_DMSO01 <- as.numeric(as.character(data_DMSO01$Count_DMSO01))
data_DMSO01 <- aggregate(data_DMSO01$Count_DMSO01, by = list(data_DMSO01$Barcode), FUN=sum)
colnames(data_DMSO01) <- c("Barcode","Count_DMSO01")

data_DMSO02 <- mydata[,c(54,53)]
colnames(data_DMSO02) <- c("Barcode","Count_DMSO02")
data_DMSO02 <- merge(x = barcodes, y = data_DMSO02[ , c("Barcode", "Count_DMSO02")], by = "Barcode", all.x = TRUE)
data_DMSO02$Barcode <- rc(data_DMSO02$Barcode)
data_DMSO02 <- data_DMSO02[order(data_DMSO02$Barcode),]
data_DMSO02$Count_DMSO02 <- as.numeric(as.character(data_DMSO02$Count_DMSO02))
data_DMSO02 <- aggregate(data_DMSO02$Count_DMSO02, by = list(data_DMSO02$Barcode), FUN=sum)
colnames(data_DMSO02) <- c("Barcode","Count_DMSO02")

data_6821 <- mydata[,c(57,56)]
colnames(data_6821) <- c("Barcode","Count_6821")
data_6821 <- merge(x = barcodes, y = data_6821[ , c("Barcode", "Count_6821")], by = "Barcode", all.x = TRUE)
data_6821$Barcode <- rc(data_6821$Barcode)
data_6821 <- data_6821[order(data_6821$Barcode),]
data_6821$Count_6821 <- as.numeric(as.character(data_6821$Count_6821))
data_6821 <- aggregate(data_6821$Count_6821, by = list(data_6821$Barcode), FUN=sum)
colnames(data_6821) <- c("Barcode","Count_6821")

data_6822 <- mydata[,c(60,59)]
colnames(data_6822) <- c("Barcode","Count_6822")
data_6822 <- merge(x = barcodes, y = data_6822[ , c("Barcode", "Count_6822")], by = "Barcode", all.x = TRUE)
data_6822$Barcode <- rc(data_6822$Barcode)
data_6822 <- data_6822[order(data_6822$Barcode),]
data_6822$Count_6822 <- as.numeric(as.character(data_6822$Count_6822))
data_6822 <- aggregate(data_6822$Count_6822, by = list(data_6822$Barcode), FUN=sum)
colnames(data_6822) <- c("Barcode","Count_6822")

data_6841 <- mydata[,c(63,62)]
colnames(data_6841) <- c("Barcode","Count_6841")
data_6841 <- merge(x = barcodes, y = data_6841[ , c("Barcode", "Count_6841")], by = "Barcode", all.x = TRUE)
data_6841$Barcode <- rc(data_6841$Barcode)
data_6841 <- data_6841[order(data_6841$Barcode),]
data_6841$Count_6841 <- as.numeric(as.character(data_6841$Count_6841))
data_6841 <- aggregate(data_6841$Count_6841, by = list(data_6841$Barcode), FUN=sum)
colnames(data_6841) <- c("Barcode","Count_6841")

data_6842 <- mydata[,c(66,65)]
colnames(data_6842) <- c("Barcode","Count_6842")
colnames(data_6842) <- c("Barcode","Count_6842")
data_6842 <- merge(x = barcodes, y = data_6842[ , c("Barcode", "Count_6842")], by = "Barcode", all.x = TRUE)
data_6842$Barcode <- rc(data_6842$Barcode)
data_6842 <- data_6842[order(data_6842$Barcode),]
data_6842$Count_6842 <- as.numeric(as.character(data_6842$Count_6842))
data_6842 <- aggregate(data_6842$Count_6842, by = list(data_6842$Barcode), FUN=sum)
colnames(data_6842) <- c("Barcode","Count_6842")

data_6881 <- mydata[,c(69,68)]
colnames(data_6881) <- c("Barcode","Count_6881")
data_6881 <- merge(x = barcodes, y = data_6881[ , c("Barcode", "Count_6881")], by = "Barcode", all.x = TRUE)
data_6881$Barcode <- rc(data_6881$Barcode)
data_6881 <- data_6881[order(data_6881$Barcode),]

data_6882 <- mydata[,c(72,71)]
colnames(data_6882) <- c("Barcode","Count_6882")
data_6882 <- merge(x = barcodes, y = data_6882[ , c("Barcode", "Count_6882")], by = "Barcode", all.x = TRUE)
data_6882$Barcode <- rc(data_6882$Barcode)
data_6882 <- data_6882[order(data_6882$Barcode),]
```

```{r}
counts_table <- data_2021[,c(1,2)]
counts_table$Count_2022 <- data_2022$Count_2022
counts_table$Count_2041 <- data_2041$Count_2041
counts_table$Count_2042 <- data_2042$Count_2042
counts_table$Count_2081 <- data_2081$Count_2081
counts_table$Count_2082 <- data_2082$Count_2082
counts_table$Count_6811 <- data_6811$Count_6811
counts_table$Count_6812 <- data_6812$Count_6812
counts_table$Count_DMSO11 <- data_DMSO11$Count_DMSO11
counts_table$Count_DMSO12 <- data_DMSO12$Count_DMSO12
counts_table$Count_DMSO21 <- data_DMSO21$Count_DMSO21
counts_table$Count_DMSO22 <- data_DMSO22$Count_DMSO22
counts_table$Count_DMSO41 <- data_DMSO41$Count_DMSO41
counts_table$Count_DMSO42 <- data_DMSO42$Count_DMSO42
counts_table$Count_DMSO81 <- data_DMSO81$Count_DMSO81
counts_table$Count_DMSO82 <- data_DMSO82$Count_DMSO82
counts_table$Count_DMSO01 <- data_DMSO01$Count_DMSO01
counts_table$Count_DMSO02 <- data_DMSO02$Count_DMSO02
counts_table$Count_6821 <- data_6821$Count_6821
counts_table$Count_6822 <- data_6822$Count_6822
counts_table$Count_6841 <- data_6841$Count_6841
counts_table$Count_6842 <- data_6842$Count_6842
counts_table$Count_6881 <- data_6881$Count_6881
counts_table$Count_6882 <- data_6882$Count_6882
#counts_table1 <- counts_table[complete.cases(counts_table), ]  ##Remove any remaining rows with NA
counts_table1 <- as.matrix(counts_table)  
rownames(counts_table1) <- counts_table1[,1]
counts_table1 <- counts_table1[,-1]
#counts_table1 <- as.matrix(sapply(counts_table1, as.numeric))  
counts_table1[is.na(counts_table1)] <- 0
storage.mode(counts_table1) <- "numeric"
```

```{r}
normalised_counts <- counts_table1
overall_median <- median(normalised_counts)
overal_mean <- mean(normalised_counts)
median_sample <- apply(normalised_counts, 2, median)
normalised_counts1 <- as.data.frame(normalised_counts)
normalised_counts1[,1] <- ( normalised_counts1[,1] / median_sample[1] ) * overall_median
normalised_counts1[,2] <- ( normalised_counts1[,2] / median_sample[2] ) * overall_median
normalised_counts1[,3] <- ( normalised_counts1[,3] / median_sample[3] ) * overall_median
normalised_counts1[,4] <- ( normalised_counts1[,4] / median_sample[4] ) * overall_median
normalised_counts1[,5] <- ( normalised_counts1[,5] / median_sample[5] ) * overall_median
normalised_counts1[,6] <- ( normalised_counts1[,6] / median_sample[6] ) * overall_median
normalised_counts1[,7] <- ( normalised_counts1[,7] / median_sample[7] ) * overall_median
normalised_counts1[,8] <- ( normalised_counts1[,8] / median_sample[8] ) * overall_median
normalised_counts1[,9] <- ( normalised_counts1[,9] / median_sample[9] ) * overall_median
normalised_counts1[,10] <- ( normalised_counts1[,10] / median_sample[10] ) * overall_median
normalised_counts1[,11] <- ( normalised_counts1[,11] / median_sample[11] ) * overall_median
normalised_counts1[,12] <- ( normalised_counts1[,12] / median_sample[12] ) * overall_median
normalised_counts1[,13] <- ( normalised_counts1[,13] / median_sample[13] ) * overall_median
normalised_counts1[,14] <- ( normalised_counts1[,14] / median_sample[14] ) * overall_median
normalised_counts1[,15] <- ( normalised_counts1[,15] / median_sample[15] ) * overall_median
normalised_counts1[,16] <- ( normalised_counts1[,16] / median_sample[16] ) * overall_median
normalised_counts1[,17] <- ( normalised_counts1[,17] / median_sample[17] ) * overall_median
normalised_counts1[,18] <- ( normalised_counts1[,18] / median_sample[18] ) * overall_median
normalised_counts1[,19] <- ( normalised_counts1[,19] / median_sample[19] ) * overall_median
normalised_counts1[,20] <- ( normalised_counts1[,20] / median_sample[20] ) * overall_median
normalised_counts1[,21] <- ( normalised_counts1[,21] / median_sample[21] ) * overall_median
normalised_counts1[,22] <- ( normalised_counts1[,22] / median_sample[22] ) * overall_median
normalised_counts1[,23] <- ( normalised_counts1[,23] / median_sample[23] ) * overall_median
normalised_counts1[,24] <- ( normalised_counts1[,24] / median_sample[24] ) * overall_median
normalised_counts1 <- normalised_counts1[,c(1,2,3,4,5,6,17,18,9,10,11,12,13,14,15,16,7,8,19,20,21,22,23,24)]
```

```{r}
normalised_countsDMSO <- normalised_counts1
normalised_countsDMSO$Baseline_Average <- (normalised_counts1$Count_DMSO01 + normalised_counts1$Count_DMSO02) / 2
normalised_countsDMSO[,1] <- normalised_countsDMSO[,1]/normalised_countsDMSO[,25]
normalised_countsDMSO[,2] <- normalised_countsDMSO[,2]/normalised_countsDMSO[,25]
normalised_countsDMSO[,3] <- normalised_countsDMSO[,3]/normalised_countsDMSO[,25]
normalised_countsDMSO[,4] <- normalised_countsDMSO[,4]/normalised_countsDMSO[,25]
normalised_countsDMSO[,5] <- normalised_countsDMSO[,5]/normalised_countsDMSO[,25]
normalised_countsDMSO[,6] <- normalised_countsDMSO[,6]/normalised_countsDMSO[,25]
normalised_countsDMSO[,7] <- normalised_countsDMSO[,7]/normalised_countsDMSO[,25]
normalised_countsDMSO[,8] <- normalised_countsDMSO[,8]/normalised_countsDMSO[,25]
normalised_countsDMSO[,9] <- normalised_countsDMSO[,9]/normalised_countsDMSO[,25]
normalised_countsDMSO[,10] <- normalised_countsDMSO[,10]/normalised_countsDMSO[,25]
normalised_countsDMSO[,11] <- normalised_countsDMSO[,11]/normalised_countsDMSO[,25]
normalised_countsDMSO[,12] <- normalised_countsDMSO[,12]/normalised_countsDMSO[,25]
normalised_countsDMSO[,13] <- normalised_countsDMSO[,13]/normalised_countsDMSO[,25]
normalised_countsDMSO[,14] <- normalised_countsDMSO[,14]/normalised_countsDMSO[,25]
normalised_countsDMSO[,15] <- normalised_countsDMSO[,15]/normalised_countsDMSO[,25]
normalised_countsDMSO[,16] <- normalised_countsDMSO[,16]/normalised_countsDMSO[,25]
normalised_countsDMSO[,19] <- normalised_countsDMSO[,19]/normalised_countsDMSO[,25]
normalised_countsDMSO[,20] <- normalised_countsDMSO[,20]/normalised_countsDMSO[,25]
normalised_countsDMSO[,21] <- normalised_countsDMSO[,21]/normalised_countsDMSO[,25]
normalised_countsDMSO[,22] <- normalised_countsDMSO[,22]/normalised_countsDMSO[,25]
normalised_countsDMSO[,23] <- normalised_countsDMSO[,23]/normalised_countsDMSO[,25]
normalised_countsDMSO[,24] <- normalised_countsDMSO[,24]/normalised_countsDMSO[,25]
normalised_countsDMSO[,17] <- normalised_countsDMSO[,17]/normalised_countsDMSO[,25]
normalised_countsDMSO[,18] <- normalised_countsDMSO[,18]/normalised_countsDMSO[,25]
normalised_countsDMSO[,25] <- normalised_countsDMSO[,25]/normalised_countsDMSO[,25]
normalised_countsDMSO <- normalised_countsDMSO[,c(25,1,2,3,4,5,6,17,18,9,10,11,12,13,14,15,16,7,8,19,20,21,22,23,24)]
```



```{r}
library(DESeq2)
library("BiocParallel")
```

```{r}
samplesheet <- read.csv(file="/home/steve/Starprom/samplesheet.csv",sep=",", row.names=1)
#rownames(counts_table1) <- counts_table1[,1]
counts <- counts_table1[, rownames(samplesheet)]
```

##Set up single and multi-factoral analysis
```{r}
analysis <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samplesheet,
                              design = ~ Drug + Time) #From sample sheet
keep <- rowSums(counts(analysis)) >= 10 ##Filter the Data
analysis <- analysis[keep,]
analysis_multi <- analysis
analysis_multi$group <- factor(paste0(analysis$Drug, analysis$Time)) ##Create multi-factoral group
group <- design(analysis_multi) <- ~ group
```

##Carry out differential analysis combining conditions
```{r include=FALSE}
analysis <- DESeq(analysis,parallel=TRUE, BPPARAM=MulticoreParam(8))
analysis_multi <- DESeq(analysis_multi,parallel=TRUE, BPPARAM=MulticoreParam(8))
```

##Single Factor Drug contrasts
```{r}
resultsNames <- resultsNames(analysis)
Lino_dmso <- results(analysis, contrast=c("Drug", "DMSO","Linopirdine"))
Diclazuril_dmso <- results(analysis, contrast=c("Drug", "DMSO","Diclazuril"))
plotMA(Lino_dmso, ylim=c(-2,2))
plotMA(Diclazuril_dmso, ylim=c(-2,2))
```

##Multi Factor Time vs Drug contrasts
```{r}
resultsNames_multi <- resultsNames(analysis_multi)
Diclazuril_time1 <- results(analysis_multi, contrast=c("group", "DiclazurilOne","DMSOOne"))
Diclazuril_time2 <- results(analysis_multi, contrast=c("group", "DiclazurilTwo","DMSOTwo"))
Diclazuril_time4 <- results(analysis_multi, contrast=c("group", "DiclazurilFour","DMSOFour"))
Diclazuril_time8 <- results(analysis_multi, contrast=c("group", "DiclazurilEight","DMSOEight"))
Lino_time2 <- results(analysis_multi, contrast=c("group", "LinopirdineTwo","DMSOTwo"))
Lino_time4 <- results(analysis_multi, contrast=c("group", "LinopirdineFour","DMSOFour"))
Lino_time8 <- results(analysis_multi, contrast=c("group", "LinopirdineEight","DMSOEight"))

plotMA(Diclazuril_time1, ylim=c(-2,2))
plotMA(Diclazuril_time2, ylim=c(-2,2))
plotMA(Diclazuril_time4, ylim=c(-2,2))
plotMA(Diclazuril_time8, ylim=c(-2,2))
plotMA(Lino_time2, ylim=c(-2,2))
plotMA(Lino_time4, ylim=c(-2,2))
plotMA(Lino_time8, ylim=c(-2,2))
```

```{r}
getmotifs <-  Lino_time2
getmotifs <- getmotifs[order(getmotifs$log2FoldChange),]
getmotifs <- subset(getmotifs, log2FoldChange > 2)
getmotifs <- getmotifs[order(getmotifs$pvalue),]
getmotifs <- subset(getmotifs, pvalue < 0.5 )
getmotifs <- head(getmotifs,20)
```

```{r}

write.table(getmotifs,"/home/steve/Starprom/lino_2_up_2lf.csv",sep=",",row.names = TRUE,col.names = TRUE)
#write.table(counts_table1,"/home/steve/Starprom/counts_table.csv",sep=",",row.names = TRUE,col.names = TRUE)
```

```{r}

write.table(normalised_counts1,"/home/steve/Starprom/normalised_counts.csv",sep="\t",row.names = TRUE,col.names = TRUE)
#write.table(normalised_countsDMSO,"/home/steve/Starprom/DMSO_normalised_counts_tab.csv",sep="\t",row.names = TRUE,col.names = TRUE)
```


```{r}
counts_table2 <- as.data.frame(counts_table1)
counts_table2$Barcode <- rownames(counts_table1)
mymotifis <- "CCTGCCACCCGCAGCTCCGG"
mymotifisrc <- rc(mymotifis)
mymotif <- mydata[(which(mydata == mymotifisrc, arr.ind = TRUE)),]
mymotif <- mydata[4062,]
mymotif2 <- counts_table2[(which(counts_table2 == mymotifis, arr.ind = TRUE)),]
mymotif3 <- counts_table[(which(counts_table == mymotifis, arr.ind = TRUE)),]

mymotif4 <- data_6842[(which(data_6842 == mymotifis, arr.ind = TRUE)),]
```

```{r}
getmotifs2 <- as.data.frame(rownames(getmotifs))
samplenames1 <- colnames(counts_table1)
colnames(getmotifs2) <- "Barcode"
motifcounts <- merge(x = getmotifs2, y = counts_table2, by = "Barcode", all.x = TRUE)
```

```{r}
ntd <- normTransform(analysis_multi)
rld <- rlog(analysis_multi, blind=FALSE)
```

##Principal component plot of the samples with Ggplot
```{r}
library(ggplot2)
pcaData <- plotPCA(rld, intgroup="DiclazurilOne", returnData=FALSE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Drug, shape=Time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
```


```{r}
library("pheatmap")
select <- order(rowMeans(counts(analysis_multi,normalized=TRUE)),
               decreasing=TRUE)[1:20]
df <- as.data.frame(colData(analysis_multi)[,c("Drug","Time")])
ntd <- normTransform(analysis_multi)
```

##Construct heatmap for the log2(n+1) data
```{r}
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
```


```{r}
rld <- rlog(analysis_multi, blind=FALSE)
#pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
```

##Heatmap of the sample-to-sample distances
```{r}
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Light_Pulse, rld$CT, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

```{r}
logcount2021 <- log10(counts_table1$Count_2021)
```


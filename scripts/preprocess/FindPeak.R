library(GenomicRanges)
args = commandArgs(T)

rootdir = args[1]
dataset = args[2]
refergenome = args[3]

print(args)

setwd(rootdir)
    
if (refergenome=="mm9")
    {
    ###load mm9 
    
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    mm9 = TxDb.Mmusculus.UCSC.mm9.knownGene
    promoter<-promoters(mm9)
    proms.df <- data.frame(promoter@seqnames, promoter@ranges@start, promoter@ranges@start + promoter@ranges@width, promoter@elementMetadata$tx_name)
    }else if (refergenome=="hg38")
    {

    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    hg38 = TxDb.Hsapiens.UCSC.hg38.knownGene
    promoter<-promoters(hg38)
    proms.df <- data.frame(promoter@seqnames, promoter@ranges@start, promoter@ranges@start + promoter@ranges@width, promoter@elementMetadata$tx_name)
    }else if (refergenome=="hg19")
    {

    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    hg19 = TxDb.Hsapiens.UCSC.hg19.knownGene
    promoter<-promoters(hg19)
    proms.df <- data.frame(promoter@seqnames, promoter@ranges@start, promoter@ranges@start + promoter@ranges@width, promoter@elementMetadata$tx_name)
    }else
    proms.df <- read.table("./mm10.refSeq_promoter.bed")


peaks.df <- read.table(paste(rootdir,dataset,"_peaks.narrowPeak",sep = ""))

# remove top 5% peaks
cutoff <- quantile((peaks.df$V5), probs = 0.95)
peaks.df <- peaks.df[which(peaks.df$V5 < cutoff),]
peaks.gr <- GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]))
proms.gr <- GRanges(proms.df[,1], IRanges(proms.df[,2], proms.df[,3]))

peaks.sel.gr <- peaks.gr[-queryHits(findOverlaps(peaks.gr, proms.gr))]
peaks.sel.ex.gr <- resize(reduce(resize(peaks.sel.gr, 1000, 
                          fix="center")), 1000, fix="center")

peaks.sel.ex.df <- as.data.frame(peaks.sel.ex.gr)[,1:3]
write.table(peaks.sel.ex.df, file = paste(rootdir,dataset,".ygi",sep = ""), 
			  append = FALSE, quote = FALSE, sep = "\t", 
			  eol = "\n", na = "NA", dec = ".", 
			  row.names = FALSE, col.names = FALSE, 
			  qmethod = c("escape", "double"),
			  fileEncoding = "")
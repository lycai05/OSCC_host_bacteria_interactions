# This is the working directory absolute path where you have the RDS files from Methylation pipeline
workdir1 = ''

# input sample files
input1 = list.files(path = workdir1, pattern = "HNmethy.*methylRaw.RDS", recursive = TRUE)
inputfiles = paste0(workdir1, input1)

# load RDS files
rds.objs = lapply(inputfiles, readRDS)
mrl <- new("methylRawList", rds.objs)
mrl@treatment <- rep(c(0,1), 44)
mrl <- filterByCoverage(mrl,lo.count=10)

# overall methylation status
par(mfrow = c(5, 10), mar = c(2, 0, 0.25, 0.25))
pdf("overallMethylStats.pdf", width = 12, height = 12)
lapply(1:length(mrl), function(x){getMethylationStats(mrl[[x]],plot=TRUE,both.strands=FALSE)})
dev.off()

# overall coverage status
par(mfrow = c(5, 10), mar = c(2, 0, 0.25, 0.25))
pdf("overallCovStats.pdf", width = 12, height = 12)
lapply(1:length(mrl), function(x){getCoverageStats(mrl[[x]],plot=TRUE,both.strands=FALSE)})
dev.off()

# Targeted CpG sites annotation
ref_peakAnno <- annotatePeak(all_CpG, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("refCpGAnno.pdf", width = 12, height = 12)
plotAnnoPie(ref_peakAnno)
dev.off()

# merge methylation status from all samples
meth=unite(mrl, destrand=FALSE, echo=FALSE, mc.cores = 16, min.per.group = 22L)

#################################################
##### Hyper DMC
#################################################

# get hyper DMC
# batchInfo=data.frame(batch = c(rep("batch1", 8), rep("batch2", 10)))
# batchAs=assocComp(mBase=meth,batchInfo)

myDiff_chisq=calculateDiffMeth(meth, mc.cores = 16, overdispersion = "MN", test = "Chisq")
myDiff = myDiff_chisq
myDiff25p.hyper=getMethylDiff(myDiff,difference=20,qvalue=0.05,type="hyper")

# get hypo DMC
myDiff25p.hypo=getMethylDiff(myDiff,difference=20,qvalue=0.05,type="hypo")

# annotate hyper DMC
#myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

myDiff25p.hyper.gr = makeGRangesFromDataFrame(myDiff25p.hyper, keep.extra.columns = F)
peakAnno <- annotatePeak(myDiff25p.hyper.gr, tssRegion=c(-3000, 3000),
                          TxDb=txdb, annoDb="org.Hs.eg.db")
hyper_peakAnno_df <- as.data.frame(peakAnno)
# annotate hypo DMC  meth=unite(mrl, destrand=FALSE, echo=FALSE, mc.cores = 16, min.per.group = 25)

#################################################
##### Hypo DMC
#################################################

# get hyper DMC
# batchInfo=data.frame(batch = c(rep("batch1", 8), rep("batch2", 10)))
# batchAs=assocComp(mBase=meth,batchInfo)

# myDiff_chisq=calculateDiffMeth(meth, mc.cores = 16, overdispersion = "MN", test = "Chisq")
myDiff25p.hypo.gr = makeGRangesFromDataFrame(myDiff25p.hypo, keep.extra.columns = F)
hypo_peakAnno <- annotatePeak(myDiff25p.hypo.gr, tssRegion=c(-3000, 3000),
                              TxDb=txdb, annoDb="org.Hs.eg.db")
hypo_peakAnno_df <- as.data.frame(hypo_peakAnno)


#################################################
##### Hyper DMR
#################################################

# identify hyper DMR
tiles=tileMethylCounts(mrl,win.size=1000,step.size=1000, mc.cores=16)
tiles_meth=unite(tiles, destrand=FALSE, mc.core = 16)

# filter metadata based on order of inputfiles
# metadata_f <- subset(metadata, !grepl("25|26|35|36|22|50|2052", methyl))
# metadata_order <- tiles_meth@sample.ids
# metadata_order <- unlist(strsplit(metadata_order, "_"))
# metadata_order <- metadata_order[grepl("HNm", metadata_order)]

# # manually change the name of the first four samples
# metadata_order[81:82] <- "HNmethy001"
# metadata_order[83:84] <- "HNmethy002"
# metadata_order[85:86] <- "HNmethy003"
# metadata_order[87:88] <- "HNmethy004"

# # prepare covariates 
# metadata_covariates <- metadata[match(metadata_order, metadata$methyl),]
# metadata_covariates <- metadata_covariates[, c("Gender", "age_c", "smoking", "Alcohol")]

# calculate DMR
tiles_myDiff=calculateDiffMeth(tiles_meth,mc.cores=16, covariates = metadata_covariates)
saveRDS(tiles_myDiff, "tiles_myDiff.rds")
# identify hyper DMR, with 25% difference
tiles_myDiff25p.hyper=getMethylDiff(tiles_myDiff,difference=20,qvalue=0.05,type="hyper")
#tiles_myDiff25p=getMethylDiff(tiles_myDiff,difference=25,qvalue=0.01)

# annotate hyper DMR, for plotting pie chart
tiles_myDiff25p.hyper.gr = makeGRangesFromDataFrame(tiles_myDiff25p.hyper, keep.extra.columns = F)
tiles_peakAnno <- annotatePeak(tiles_myDiff25p.hyper.gr, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")

# get hyper DMR df
tiles_peakAnno_df <- as.data.frame(tiles_peakAnno)
tiles_peakAnno_df <- merge(tiles_peakAnno_df, getData(tiles_myDiff25p.hyper), by.y = c("chr", "start", "end"), by.x = c("seqnames", "start", "end"))
tiles_peakAnno_df <- tiles_peakAnno_df[order(tiles_peakAnno_df$qvalue), ]


tile_hyper_promoters <- subset(tiles_peakAnno_df, grepl('Promoter', annotation))
tile_hyper_promoters <- tile_hyper_promoters[order(-tile_hyper_promoters[["qvalue"]], tile_hyper_promoters[["meth.diff"]], decreasing = T), ]


tile_hyper_promoters_merged <- merge(tile_hyper_promoters, getData(tiles_meth), by.x = c("seqnames", "start", "end"), by.y = c("chr", "start", "end"), all.x = T)

sink("tile_hyper_promoters.txt")
cat(tile_hyper_promoters$SYMBOL, sep = "\n")
sink()

write.table(tile_hyper_promoters_merged, "tile_hyper_promoters_merged.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#################################################
##### Hypo DMR
#################################################
# identify hypo DMR, with 25% dfference
tiles_myDiff25p.hypo=getMethylDiff(tiles_myDiff,difference=20,qvalue=0.05,type="hypo")

# annotate hypo DMR, for plotting pie chart
tiles_myDiff25p.hypo.gr = makeGRangesFromDataFrame(tiles_myDiff25p.hypo, keep.extra.columns = F)
tiles_hypo_peakAnno <- annotatePeak(tiles_myDiff25p.hypo.gr, tssRegion=c(-3000, 3000),
                                    TxDb=txdb, annoDb="org.Hs.eg.db")

# get hypo DMR df
tiles_hypo_peakAnno_df <- as.data.frame(tiles_hypo_peakAnno)
tiles_hypo_peakAnno_df <- merge(tiles_hypo_peakAnno_df, getData(tiles_myDiff25p.hypo), by.y = c("chr", "start", "end"), by.x = c("seqnames", "start", "end"))
tiles_hypo_peakAnno_df <- tiles_hypo_peakAnno_df[order(tiles_hypo_peakAnno_df$qvalue), ]

tile_hypo_promoters <- subset(tiles_hypo_peakAnno_df, grepl('Promoter', annotation))
tile_hypo_promoters <- tile_hypo_promoters[order(tile_hypo_promoters[["qvalue"]], tile_hypo_promoters[["meth.diff"]], decreasing = F), ]
tile_hypo_promoters_merged <- merge(tile_hypo_promoters, getData(tiles_meth), by.x = c("seqnames", "start", "end"), by.y = c("chr", "start", "end"), all.x = T)

sink("tile_hypo_promoters.txt")
cat(tile_hypo_promoters$SYMBOL, sep = "\n")
sink()

write.table(tile_hypo_promoters_merged, "tile_hypo_promoters_merged.txt", sep = "\t", row.names = F, col.names = T, quote = F)
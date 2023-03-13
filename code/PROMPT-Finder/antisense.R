library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(rstatix)
library(ggpubr) 
library(eulerr)
normalize_length = "TRUE"
select_length = "TRUE"
# the cutoff of high expression level in ctrl
args = commandArgs(T)
cut <- 20
site_cut <- 500
treated = args[1] #"IAA"
species = args[2] # "hg38"

######### functions ############
volcano <- function(df, outname, ylimit=FALSE, xlimit=FALSE, fc=0){
  pval_threshold <- 0.05
  logfc_threshold <- fc
  color_factor <- as.factor(ifelse(df$log2FoldChange > logfc_threshold & 
                                     df$padj <= pval_threshold, "a",
                                   ifelse(df$log2FoldChange < -logfc_threshold & 
                                            df$padj <= pval_threshold, "b", "c")))
  color_list = c("red", "blue","gray")
  #color_list = c("red","gray")
  #print(color_list)
  significance_up <-df %>% filter(log2FoldChange >=logfc_threshold & padj <= pval_threshold) %>% nrow()
  significance_down <-df %>% filter(log2FoldChange <= -logfc_threshold & padj <= pval_threshold) %>% nrow()
  volcano = ggplot(data=df, 
                   aes(x=log2FoldChange, y=-log10(padj), 
                       colour=color_factor)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_classic() +
    scale_color_manual(values = color_list)+
    xlab("log2 fold change") + ylab("-log10 FDR") +
    labs(title = sub("(.*).pdf", "\\1", outname), 
         subtitle = paste0("n=", nrow(df), "   significance_up=", significance_up, "   significance_down=", significance_down)) +
    theme(legend.position = 'none') 
  #volcano = volcano + geom_text_repel(aes(label=ifelse(padj <= 0.05 & abs(log2FoldChange) >= logfc_threshold,name2, '')), max.overlaps = 15)
  if (min(df$log2FoldChange, na.rm = TRUE) >= 0) {
    volcano = volcano + xlim(-1, max(df$log2FoldChange, na.rm = TRUE)) + geom_point(color="red")
  }
  if (ylimit != FALSE) {
    volcano = volcano + coord_cartesian(ylim=c(0, ylimit))
  }
  if (xlimit[1] != FALSE) {
    volcano = volcano + xlim(xlimit)
  }
  ggsave(outname, width = 1800,height = 1800, dpi = 300, units = "px", plot = volcano)
}

boxplot_n <- function(counts_temp, n){
  counts_temp_wider = counts_temp %>% pivot_wider(names_from = samples, values_from = counts) %>% dplyr::select(-ID)
  error_site = sapply(counts_temp_wider,function(x) boxplot.stats(x)$stats[5])
  ymax = max(error_site)
  stat.test <- counts_temp  %>% wilcox_test(counts ~ samples) %>% adjust_pvalue() %>% add_significance("p.adj") %>% add_x_position(x = "samples")
  stat.test$y.position = seq(ymax*1.1, ymax*1.1 + (ymax*0.1) * (nrow(stat.test) - 1), ymax*0.1)
  p <- ggplot(counts_temp, aes(x = as.factor(samples), y = counts)) + ylab("Normalized counts") + xlab("") +
    stat_boxplot(geom='errorbar', linetype=1, width=0.2, aes(color = samples))+
    geom_boxplot(width = 0.5,  outlier.colour = NA, aes(color = samples)) +
    theme_classic()+theme(legend.position = "none") + labs(subtitle = paste0("n=", nrow(counts_temp) / length(unique(counts_temp$samples))))
  p + coord_cartesian(ylim =c(0, ymax*1.3)) + stat_pvalue_manual(stat.test, label = "p.adj",tip.length = 0)
  ggsave(paste0(n, "_boxplot.pdf"), dpi = 300, width = 1300, height = 2000, units = "px")
}

boxplot_2_counts<- function(counts, mark){
  # plot with counts ----sense
  counts_temp = counts %>%  pivot_longer(cols = -ID, values_to = "counts", names_to = "samples") 
  counts_temp$samples = relevel(as.factor(counts_temp$samples),ref=control)
  boxplot_n(counts_temp, mark)
}

venn <- function(ls, name){
  a <- plot(euler(ls, shape = "ellipse"), quantities = TRUE, main = name)
  ggsave(name, plot = a, width = 5, height = 5)
}

{control = "CTRL"
  database <- read.csv('database.csv',header = T,check.names = FALSE)
  ##arrange the database, same as the range of featurecounts
  database_rep <- database %>% group_by(condition) %>% mutate(su = length(name)) %>% mutate(rep = 1:unique(su)) %>% 
    dplyr::select(-su) %>% mutate(condition_rep = paste(condition,"rep",rep,sep = "_")) %>% mutate(name = condition_rep)
  data <- read.table('featureCounts.txt',header = T, quote = '\t', skip = 1) 
  data[,7:ncol(data)] <- lapply(data[,7:ncol(data)], as.integer)
  Len <- data %>% dplyr::select(transcriptID = Geneid, Length)
  sampleNmaes <- as.factor(database_rep$condition_rep)
  countData <-as.matrix(data[,7:ncol(data)])
  colnames(countData) <- sampleNmaes
  rownames(countData) <- data$Geneid
  rownames(database) <- sampleNmaes
  database$condition <- as.factor(database$condition)
  database$condition <- relevel(database$condition,ref=control)
  if (length(args) == 3) {
    spikein="T"
    # size factor = 1/scale factor
    sf_df = read.csv(args[3], sep = "\t")
    sf_df$name = str_replace_all(sf_df$sample, ".bam", "") %>% as.numeric()
    sf_df = left_join(database, sf_df, by="name")
    sizefactor = 1 / sf_df$scalingFactor
    print("sizefactor")
    print(sizefactor)
  } else{
    spikein="F"
  } 
  ##DESeq2
  dds <- DESeqDataSetFromMatrix(countData, colData = database, design = ~condition)
  mcols(dds)$basepairs <- Len[,"Length"]
  if (spikein=="T") {
    # normalized with spike-in's size factor
    dds = estimateSizeFactors(dds)
    sizeFactors(dds) = sizefactor
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)
  }else{
    dds <- DESeq(dds)
  }
  vsd = vst(dds, blind=FALSE)
  vst.dds <- vst(dds)
  vst <- assay(vst.dds)
}
sf = sizeFactors(dds) %>% as.data.frame() %>% rownames_to_column("sample")
colnames(sf) = c("sample", "scalingFactor")
sf$sample = str_replace(sf$sample, "_rep_.", "")
sf$scalingFactor = 1 / sf$scalingFactor 
write.table(sf, file = "scalingfactor.txt", sep = "\t", quote = FALSE, row.names = FALSE)

fpkm_ = fpkm(dds) 
gene_fpkm = fpkm_ %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(!startsWith(ID, "RT_"))
ua_fpkm = fpkm_ %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_"))
# keep all ua(PROMPTs)
all_ua = ua_fpkm$ID

for (i in unique(database$condition)) {
  temp = gene_fpkm %>% dplyr::select(contains(i)) %>% rowMeans() %>% as.data.frame()
  colnames(temp) = i
  gene_fpkm = gene_fpkm %>% cbind(temp)  
  temp = ua_fpkm %>% dplyr::select(contains(i)) %>% rowMeans() %>% as.data.frame()
  colnames(temp) = i
  ua_fpkm = ua_fpkm %>% cbind(temp) 
}
gene_fpkm_all = gene_fpkm
gene_fpkm = gene_fpkm %>% dplyr::select(!contains("rep")) %>% filter(CTRL >1)
gene_fpkm_keep = gene_fpkm
ua_fpkm_all = ua_fpkm = ua_fpkm %>% dplyr::select(!contains("rep")) %>%
  mutate(ID = gsub("RT_(.*)", "\\1", ID))
ua_fpkm = ua_fpkm_all %>% filter(ID %in% gene_fpkm$ID)
# venn plot for select gene and PROMPTs
venn(list(all_gene=gene_fpkm_all$ID, expressed_gene=gene_fpkm_keep$ID, PROMPTs=ua_fpkm_all$ID), "selectgene.pdf")

# contrast treated and ctrl 
res = results(dds, contrast = c("condition",treated,"CTRL"))
res_gene = res %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(!startsWith(ID, "RT_"))
res_ua = res %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_"))
res_ua$ID = gsub("RT_(.*)", "\\1", res_ua$ID)
# volcanno plot for PROMPTs 
volcano(res_ua %>% filter(ID %in% gene_fpkm_keep$ID), "ua_vol.pdf", ylimit = 200, fc=1)

volcano(res_ua %>% filter(ID %in%  str_replace(all_ua, "RT_", "")) , "uaall_vol.pdf", ylimit = 200, fc=1)



##get TSS distance PN df 
{
  ##filter antisense TSS site and filter TSS site at 0-500 bp
  NM = read.table(paste0("~/reference/",species,"/RefSeq.bed"),sep = "\t") %>%  filter(!str_detect(V1,"_"))
  RT = read.table("merge.bed",sep = "\t") %>% filter(!str_detect(V1,"_")) %>% mutate(V4=gsub(".*_(N.*)","\\1",V4))
  # calculate length distribution of PROMPTs
  len = RT %>% mutate(length = V3-V2) 
  x <- ggplot(len,aes(length))+
    geom_histogram(binwidth = 50)+
    theme(legend.position = 'none') +
    xlab("distance(bp)")+ 
    xlim(c(0, 20000))+
    labs(title = "The distribution of PROMPTs length" ,subtitle = paste0("n=", nrow(len)))+
    theme_classic() + theme()
  ggsave('PROMPTslenth.pdf',width = 2400,height = 1800, dpi = 300, units = "px", plot = x)
  NM_p <- NM %>% 
    filter(V6 == "+") %>% 
    dplyr::select(V4,V2)
  NM_m = NM %>% 
    filter(V6 == "-") %>% 
    dplyr::select(V4,V3)
  RT_p <-RT %>% 
    filter(V6 == "-") %>% 
    dplyr::select(V4,V3)
  RT_m <-RT %>% 
    filter(V6 == "+") %>% 
    dplyr::select(V4,V2)
  Positive <- left_join(NM_p,RT_p,by = "V4") %>% mutate(site = V2 - V3) %>% dplyr::select(V4,site)
  Negavite <- left_join(NM_m,RT_m,by="V4") %>% mutate(site = V2-V3) %>%  dplyr::select(V4,site)
  PN <- rbind(Positive,Negavite) %>%  filter(site < site_cut & site > 0)  # %>% filter(V4 %in% name_list)
  PN_list <- PN$V4  ##the list of filtered transcript id
}

# site distribution
PN2 <- rbind(Positive,Negavite)  %>% distinct(V4, .keep_all = TRUE) %>% 
  filter(site >=-2000 &site<= 2000) %>% filter(V4 %in% str_replace(all_ua, "RT_", ""))
PN2_300<-PN2 %>% filter(site>0 & site < site_cut) %>% distinct(V4, .keep_all = TRUE)
list_300 <- PN2_300$V4
write.table(PN2_300$V4,"list300.txt",quote = FALSE, row.names = FALSE,col.names = FALSE)
x <- ggplot(PN2,aes(-site))+
  geom_histogram(binwidth = 10)+
  theme(legend.position = 'none') +
  xlab("distance(bp)")+ 
  xlim(c(-2000, 2000))+
  labs(title = "The distribution of antisense TSS sits" ,subtitle = paste0("n=", nrow(PN2), "   0-", site_cut, "(bp)=", nrow(PN2_300)))+
  theme_classic() + theme()
x
ggsave('RT_TSS_site.pdf',width = 2400,height = 1800, dpi = 300, units = "px", plot = x)

# boxplot for fpkm
res_uaup = res_ua %>% filter(log2FoldChange >1, padj <0.05)
write.csv(res_uaup, "uaup.csv", quote = FALSE, row.names = FALSE)
ua_fpkm = ua_fpkm %>% filter(ID %in% res_uaup$ID)
gene_fpkm = gene_fpkm %>% filter(ID %in% ua_fpkm$ID)
boxplot_2_counts(ua_fpkm, "PROMPT_fpkm")
boxplot_2_counts(gene_fpkm, "gene_fpkm")
gene_bed = NM %>% filter(V4 %in% gene_fpkm$ID)
write.table(gene_bed, "gene_bed.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# generate normalized counts with size facotor and transcripts length

count_ = counts(dds, normalized=TRUE) %>% merge(Len %>% column_to_rownames("transcriptID"), by=0) %>% column_to_rownames("Row.names")
gene_count = count_ %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(!startsWith(ID, "RT_"))
ua_count = count_ %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(startsWith(ID, "RT_"))

for (i in unique(database$condition)) {
  temp = gene_count %>% dplyr::select(contains(i)) %>% rowMeans() %>% as.data.frame()
  colnames(temp) = i
  gene_count = gene_count %>% cbind(temp)  
  temp = ua_count %>% dplyr::select(contains(i)) %>% rowMeans() %>% as.data.frame()
  colnames(temp) = i
  ua_count = ua_count %>% cbind(temp) 
}
gene_count_temp = gene_count %>% select(!contains("rep_")) %>% column_to_rownames("ID") %>% select(-Length)
gene_count = (gene_count_temp * 1000 / gene_count$Length) %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(ID %in% gene_fpkm$ID)
ua_count_temp = ua_count %>% select(!contains("rep_")) %>% column_to_rownames("ID") %>% select(-Length)
rownames(ua_count_temp) = gsub("RT_(.*)", "\\1", rownames(ua_count_temp))
ua_count = (ua_count_temp * 1000 / ua_count$Length) %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(ID %in% ua_fpkm$ID)
boxplot_2_counts(ua_count, "PROMPT_nc")
boxplot_2_counts(gene_count, "gene_nc")



# quantile expression in ctrl 
gene_fpkm$quantile = ntile(gene_fpkm$CTRL, 4)
write.csv(gene_fpkm, "gene_fpkm.csv")

# U1 score for low and high expression score
u1_score = read.csv("u1_gene_score.txt", sep="\t") %>% select(ID = name, scores)

gene_fpkm = gene_fpkm %>% left_join(u1_score, by = "ID")
boxplot_n(gene_fpkm %>% select(ID, samples=quantile, counts=scores), "score")

res_gene  = res_gene %>% left_join(u1_score, by="ID")

# sort gene with up down regulate after IAA treated
res_gene$change = ifelse(res_gene$log2FoldChange < 1, "down", "up")
# boxplot for fc & score
boxplot_n(res_gene %>% select(ID, samples=change, counts=scores) %>% na.omit(), "changescore")
# point and correlated for fc & score
ggplot(res_gene, aes(x=log2FoldChange, y=scores)) + geom_point(alpha=0.8, color="gray4") + theme_bw() + 
  labs(subtitle = paste0("pearson = ", cor(res_gene$scores, res_gene$log2FoldChange, use = "complete.obs", method = "pearson"),
                         "\nn = ", nrow(res_gene)))
ggsave("changescore_vol.pdf", width = 8, height = 8)

res_gene = res_gene %>% left_join(gene_fpkm_all, by="ID")
boxplot_n(res_gene %>% select(ID, samples=change, counts=CTRL) %>% na.omit(), "changefpkm")
# point and correlated for fc & fpkm
ggplot(res_gene, aes(x=log2FoldChange, y=log2(CTRL+1))) + geom_point(alpha=0.8, color="gray4") + theme_bw() + 
  labs(subtitle = paste0("pearson = ", cor(res_gene$CTRL, res_gene$log2FoldChange, use = "complete.obs", method = "pearson"),
                         "\nn = ", nrow(res_gene)))
ggsave("changefpkm_vol.pdf", width = 8, height = 8)

# # fpkm and u1 score
# res_gene$quantile = ntile(res_gene$CTRL, 4)
# 
# boxplot_n(res_gene %>% select(ID, samples=quantile, counts=scores) %>% na.omit(), "fpkmscore")
# # point and correlated for fc & fpkm
# ggplot(res_gene, aes(x=log2(CTRL+1), y=scores)) + geom_point(alpha=0.8, color="gray4") + theme_bw() + 
#   labs(subtitle = paste0("pearson = ", cor(res_gene$CTRL, res_gene$scores, use = "complete.obs", method = "pearson"),
#                          "\nn = ", nrow(res_gene)))
# ggsave("fpkmscore_vol.pdf", width = 8, height = 8)


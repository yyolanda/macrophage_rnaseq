library(plyr)
library(dplyr)
library(data.table)
library(limma)
library(edgeR)
library(readxl)
library(stringr)
source('/data/CODE/DMR_functions.r')
setwd('/data/Don/RNA_Seq_Projects/Projects/Hiroto_RNAseq')
options(stringsAsFactors=F)


############################ GENE LEVEL ANALYSIS #########################################
# load gene counts
gene_cnt = readRDS('processed_data/gene_expectedCount_RSEM.rds')
SampleName = gsub('(BEC|BIDC)-([0-9]+).*','\\1\\2',colnames(gene_cnt))[-c(1)]
MType = gsub('(BEC|BIDC)-([0-9]+)-(DP|DN|M1|M2)-RNA.*','\\3',colnames(gene_cnt))[-c(1)]
colnames(gene_cnt) = gsub('(BEC|BIDC)-([0-9]+)-(DP|DN|M1|M2)-RNA.*','\\1_\\2_\\3',colnames(gene_cnt))
# load phenotypes
pDat = fread('pheno.csv')
RIN = read_xlsx('RIN and methods.xlsx')
RIN$SampleName=gsub(' RNA','',RIN$SampleName)
RIN$SampleName=gsub(' ','_',RIN$SampleName)
pDat = merge(pDat, RIN, by.x='FileName',by.y='SampleName', all=T)
pDat$Method=recode(pDat$Method,`RLT puls`='RLT',`PBS+10%FBS`='PBS')
anno = fread('processed_data/gtf_clean.csv')

# normalization and filtering
id = gene_cnt$gene_id
gene_cnt=gene_cnt[,-c(1)]
gene_cnt=as.data.frame(gene_cnt)
rownames(gene_cnt)=id

# Standard usage of limma/voom
design = model.matrix(~MType, data=pDat)
dge <- DGEList(counts=gene_cnt)
keep <- filterByExpr(dge,design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE);dev.off()
# re-normalize with the intra-cluster correlation
dupcor <- duplicateCorrelation(v,design,block=pDat$SampleName)
vobj = voom(dge, design, plot=TRUE, block=pDat$SampleName, correlation=dupcor$consensus);dev.off() #16452
saveRDS(vobj$E,'processed_data/gene_normfilt_RSEM.rds')

# load normalized data
eDat = readRDS('processed_data/gene_normfilt_RSEM.rds')
eDat = eDat[rowSums(eDat > 1)> 37/4,] #12641
pDat = fread('pheno.csv')
RIN = read_xlsx('RIN and methods.xlsx')
RIN$SampleName=gsub(' RNA','',RIN$SampleName)
RIN$SampleName=gsub(' ','_',RIN$SampleName)
pDat = merge(pDat, RIN, by.x='FileName',by.y='SampleName', all=T)
pDat$Method=recode(pDat$Method,`RLT puls`='RLT',`PBS+10%FBS`='PBS')
anno = fread('processed_data/gtf_clean.csv')


########################## PCA #############################
# PCA by macrophage subtype
pca <- prcomp(t(eDat), scale. = TRUE)
Group = pDat$MType
table(Group)
#DN DP M1 M2
# 9  8 10 10
p1=ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = Group,
         ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
         circle = FALSE) + theme_bw()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c("#e41a1c", "#377eb8",'#4daf4a','#984ea3'))+
  theme(axis.title=element_text(face="bold",size="14"),
		axis.text=element_text(face="bold",size="14"),
		legend.text=element_text(face="bold",size="11"),
		legend.title=element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave('results/PCA_mtype.png', p1, height = 5, width = 5)

# PCA by RIN
Group = ifelse(pDat$RIN<=8,'<=8','>8')
table(Group)
p2=ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = Group,
         ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
         circle = FALSE) + theme_bw()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c("#e41a1c", "#377eb8",'#4daf4a','#984ea3'))+
  theme(axis.title=element_text(face="bold",size="14"),
		axis.text=element_text(face="bold",size="14"),
		legend.text=element_text(face="bold",size="11"),
		legend.title=element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave('results/PCA_RIN.png', p2, height = 5, width = 5)

# PCA by smoking status
Group = pDat$SMOKING_STATUS
table(Group)
p3=ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = Group,
         ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
         circle = FALSE) + theme_bw()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c("#e41a1c", "#377eb8",'#4daf4a','#984ea3'))+
  theme(axis.title=element_text(face="bold",size="14"),
		axis.text=element_text(face="bold",size="14"),
		legend.text=element_text(face="bold",size="11"),
		legend.title=element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave('results/PCA_Smoking.png', p3, height = 5, width = 5)

# PCA by disease
Group = pDat$Disease
table(Group)
p4=ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = Group,
         ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
         circle = FALSE) + theme_bw()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c("#e41a1c", "#377eb8",'#4daf4a','#984ea3'))+
  theme(axis.title=element_text(face="bold",size="14"),
		axis.text=element_text(face="bold",size="14"),
		legend.text=element_text(face="bold",size="11"),
		legend.title=element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave('results/PCA_disease.png', p4, height = 5, width = 5)

# PCA by method
Group = pDat$Method
table(Group)
p5=ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = Group,
         ellipse = TRUE, varname.size=TRUE, pc.biplot=FALSE, var.axes=FALSE,
         circle = FALSE) + theme_bw()+
  theme(legend.position = 'bottom')+
  scale_color_manual(values = c("#e41a1c", "#377eb8",'#4daf4a','#984ea3'))+
  theme(axis.title=element_text(face="bold",size="14"),
		axis.text=element_text(face="bold",size="14"),
		legend.text=element_text(face="bold",size="11"),
		legend.title=element_blank())+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
ggsave('results/PCA_method.png', p5, height = 5, width = 5)

ggdraw() +
  draw_plot(p1, x = 0.01, y = .50, width = 0.48, height = .48) +
  draw_plot(p2, x = 0.52, y = .50, width = 0.48, height = .48) +
  draw_plot(p3, x = 0.01, y = .01, width = 0.48, height = .48) +
  draw_plot(p5, x = 0.52, y = .01, width = 0.48, height = .48) +
  draw_plot_label(label = c("A","B","C","D"), size = 15,
                  x = c(0, 0.5, 0,0.5), y = c(1,1,0.5,0.5))
ggsave('results/PCAs.png',height = 10,wid=10)

##################### Differential expression ######################
# cell-specific markers
design <- model.matrix(~0+MType+Method+SMOKING_STATUS, data=pDat)
corfit <- duplicateCorrelation(eDat,design,block=pDat$SampleName)
corfit$consensus #intra-subject correlation
fit <- lmFit(eDat,design,block=pDat$SampleName,correlation=corfit$consensus)
cm <- makeContrasts(
		MTypeDN-(MTypeDP+MTypeM1+MTypeM2)/3,
		MTypeDP-(MTypeDN+MTypeM1+MTypeM2)/3,
		MTypeM1-(MTypeDN+MTypeDP+MTypeM2)/3,
		MTypeM2-(MTypeDN+MTypeDP+MTypeM1)/3,
		levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

topt = topTable(fit2, coef=1, num=Inf)
topt_annot = topt %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
fwrite(topt_annot, 'results/DN.csv')

topt = topTable(fit2, coef=2, num=Inf)
topt_annot = topt %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
fwrite(topt_annot, 'results/DP.csv')

topt = topTable(fit2, coef=3, num=Inf)
topt_annot = topt %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
fwrite(topt_annot, 'results/M1.csv')

topt = topTable(fit2, coef=4, num=Inf)
topt_annot = topt %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
fwrite(topt_annot, 'results/M2.csv')

dT = decideTests(fit2,method="global",adjust.method="BH",p.value=0.1,lfc=0)
dT_annot = dT %>% as.data.frame() %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
summary(dT)
#       MTypeDN - (MTypeDP + MTypeM1 + MTypeM2)/3
#Down                                         691
#NotSig                                     11675
#Up                                           275
#       MTypeDP - (MTypeDN + MTypeM1 + MTypeM2)/3
#Down                                         178
#NotSig                                     12195
#Up                                           268
#       MTypeM1 - (MTypeDN + MTypeDP + MTypeM2)/3
#Down                                          18
#NotSig                                     12562
#Up                                            61
#       MTypeM2 - (MTypeDN + MTypeDP + MTypeM1)/3
#Down                                          88
#NotSig                                     12513
#Up                                            40

topt = rbind(fread('results/DN.csv') %>% mutate(comparison = 'DN'),
			fread('results/DP.csv') %>% mutate(comparison = 'DP'),
			fread('results/M1.csv') %>% mutate(comparison = 'M1'),
			fread('results/M2.csv') %>% mutate(comparison = 'M2'))
topt = topt %>% mutate(global_FDR = p.adjust(P.Value, method='BH'))
fwrite(topt,'results/full_cellmarker_results.csv')




####################### enrichment analysis ######################
library(xlsx)
dT_annot=fread('results/full_cellmarker_results.csv')
ref = dT_annot %>% dplyr::select(Gene) %>% unlist() %>% unique()

degenes = dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'DN') %>% dplyr::select(Gene) %>% unlist() %>% unique()
DN_enrich = enrich.analysis(degenes,ref)
degenes = dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'DP') %>% dplyr::select(Gene) %>% unlist() %>% unique()
DP_enrich = enrich.analysis(degenes,ref)
degenes = dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'M1') %>% dplyr::select(Gene) %>% unlist() %>% unique()
M1_enrich = enrich.analysis(degenes,ref)
degenes = dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'M2') %>% dplyr::select(Gene) %>% unlist() %>% unique()
M2_enrich = enrich.analysis(degenes,ref)

GOBP_enrich_tab = rbind(
	if(nrow(DN_enrich$GOBP@result)!=0) DN_enrich$GOBP %>% simplify(cutoff=0.7) %>% .@result %>% dplyr::select(ID, Description, p.adjust, geneID, Count) %>%
		filter(p.adjust<0.1) %>% mutate(Category = 'DN'),
	if(nrow(DP_enrich$GOBP@result)!=0) DP_enrich$GOBP %>% simplify(cutoff=0.7) %>% .@result %>% dplyr::select(ID, Description, p.adjust, geneID, Count) %>%
		filter(p.adjust<0.1) %>% mutate(Category = 'DP'),
	if(nrow(M1_enrich$GOBP@result)!=0) M1_enrich$GOBP %>% simplify(cutoff=0.7) %>% .@result %>% dplyr::select(ID, Description, p.adjust, geneID, Count) %>%
		filter(p.adjust<0.1) %>% mutate(Category = 'M1'),
	if(nrow(M2_enrich$GOBP@result)!=0) M2_enrich$GOBP %>% simplify(cutoff=0.7) %>% .@result %>% dplyr::select(ID, Description, p.adjust, geneID, Count) %>%
		filter(p.adjust<0.1) %>% mutate(Category = 'M2'))

write.xlsx(GOBP_enrich_tab, file="results/enrichment_results_cellmarker_FCgt1.xlsx", sheetName="GOBP_nonRd", append=TRUE, row.names=FALSE)


# plot top GOBP pathways in each category
library(readxl)
dat=read_xlsx('results/enrichment_results_cellmarker_FCgt1.xlsx', sheet=3)

plotdat = as.data.table(dat)[,head(.SD %>% filter(Count>2) %>% arrange(p.adjust,Count),5),by=Category]
plotdat = plotdat %>% mutate(number=1:nrow(plotdat))
setnames(plotdat,'p.adjust','FDR')
ggplot(plotdat, aes(x=Category, y=reorder(Description,number), size=Count, fill=-log10(FDR)))+
  geom_point(shape=21)+
  scale_size_area(name = 'Number of Overlapping Genes', max_size = 10)+
  scale_fill_gradient(low = "#e8d3ff", high = "#240049", space = "Lab",
                      na.value = "grey50", guide = "colourbar",  name = expression(paste('- Log'[10],' FDR'))) +
  xlab('Macrophage Subtypes')+
  ylab('GO Biological Processes')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(axis.title=element_text(face="bold",size=11),axis.text=element_text(size=11,face="bold"))
ggsave('results/enrichment_GOBP_cellmarker.png',height=7, width = 10)





######################### volcano plots #############################
# volcano plot function
volcano_plot_exprs <-
function(toptable, contrast.fld = 'betaFC', pval.fld = 'P.Value',
                         pthresh = 0.01, contrast.thresh = 0.25,
                         xlabel = 'log2 fold change',
                         ylabel = '-log10 p-value',
                         title = '',
                         label.fld = 'Gene',
                         line.label='',
                         label.thresh = 0.05, label.thresh.fld = 'adj.P.Val',
                         colors = c(neut = '#8b8b8b', nonsig='grey', low = 'blue', high = 'red')){
  require(ggplot2);require(ggrepel);require(ggrastr)

  df <- data.frame(contrast = toptable[,contrast.fld],
                   pval = toptable[,pval.fld])
  df$colorClass <- 'nonsig'
  df$colorClass[df$pval < pthresh] <- 'neut'
  df$colorClass[df$pval < pthresh & df$contrast < -contrast.thresh] <- 'Down-Regulated'
  df$colorClass[df$pval < pthresh & df$contrast > contrast.thresh] <- 'Up-Regulated'
  names(colors)[3:4] <- c('Down-Regulated', 'Up-Regulated')

  ### labeling:
  df$label <- ''
  df$hjust <- as.numeric(df$contrast < 0)
  if(!is.null(label.fld)){
    ix <- which(toptable[[label.thresh.fld]] < label.thresh& toptable[[contrast.fld]]>0)
    df$label[ix] <- toptable[[label.fld]][ix]
  }

  ### add label to p-value threshold line:
  pval.line.label <- line.label
  pval.x <- min(toptable[[contrast.fld]])

  ggplot(df, aes(x = contrast, y = -log10(pval), color = colorClass)) +
    geom_point_rast(size = 2) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(name = 'Effect Direction', values = colors) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text_repel(aes(label=label), size = 3.5) +
    theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'),
          plot.title = element_text(size = 15, face = 'bold'),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = c(0,0), legend.justification = c(0,0),
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
             hjust = 0, vjust = 0, size = rel(4), fontface = 3)
}


pcutoff=dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'DN') %>% dplyr::select(P.Value) %>% unlist() %>% max()
p=volcano_plot_exprs(dT_annot %>% filter(comparison == 'DN'),
					contrast.fld = 'logFC', pval.fld = 'P.Value',
                    pthresh = pcutoff, contrast.thresh = 0,
                    xlabel = expression(paste(log[2], ' Fold Change')),
                    ylabel = expression(paste('-',log[10], ' P')),
                    label.fld = 'Gene',
                    line.label='FDR=0.1',
                    label.thresh = 1e-4, label.thresh.fld = 'adj.P.Val')+
	theme(legend.position='none')
ggsave('results/volcano_DN.png', height =5, width=5)

pcutoff=dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'DP') %>% dplyr::select(P.Value) %>% unlist() %>% max()
p=volcano_plot_exprs(dT_annot %>% filter(comparison == 'DP'),
					contrast.fld = 'logFC', pval.fld = 'P.Value',
                    pthresh = pcutoff, contrast.thresh = 0,
                    xlabel = expression(paste(log[2], ' Fold Change')),
                    ylabel = expression(paste('-',log[10], ' P')),
                    label.fld = 'Gene',
                    line.label='FDR=0.1',
                    label.thresh = 5e-3, label.thresh.fld = 'adj.P.Val')+
	theme(legend.position='none')
ggsave('results/volcano_DP.png', height =5, width=5)

pcutoff=dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'M1') %>% dplyr::select(P.Value) %>% unlist() %>% max()
p=volcano_plot_exprs(dT_annot %>% filter(comparison == 'M1'),
					contrast.fld = 'logFC', pval.fld = 'P.Value',
                    pthresh = pcutoff, contrast.thresh = 0,
                    xlabel = expression(paste(log[2], ' Fold Change')),
                    ylabel = expression(paste('-',log[10], ' P')),
                    label.fld = 'Gene',
                    line.label='FDR=0.1',
                    label.thresh = 0.095, label.thresh.fld = 'adj.P.Val')+
	theme(legend.position='none')
ggsave('results/volcano_M1.png', height =5, width=5)

pcutoff=dT_annot %>% filter(adj.P.Val<0.1 & comparison == 'M2') %>% dplyr::select(P.Value) %>% unlist() %>% max()
p=volcano_plot_exprs(dT_annot %>% filter(comparison == 'M2'),
					contrast.fld = 'logFC', pval.fld = 'P.Value',
                    pthresh = pcutoff, contrast.thresh = 0,
                    xlabel = expression(paste(log[2], ' Fold Change')),
                    ylabel = expression(paste('-',log[10], ' P')),
                    label.fld = 'Gene',
                    line.label='FDR=0.1',
                    label.thresh = 0.05, label.thresh.fld = 'adj.P.Val')+
	theme(legend.position='none')
ggsave('results/volcano_M2.png', height =5, width=5)




#################### WGCNA ####################################
TPM = fread('processed_data/genes_RSEM_GENCODE_count.tsv')
colnames(TPM) = gsub('(BEC|BIDC)-([0-9]+)-(DP|DN|M1|M2)-RNA.*','\\1_\\2_\\3',colnames(TPM))
rn = TPM$gene_id
TPM = TPM %>% dplyr::select(-gene_id, -transcript_id.s.) %>% as.data.frame()
rownames(TPM)=rn
saveRDS(TPM, 'processed_data/gene_TPM_RSEM.rds')

TPM = readRDS('processed_data/gene_TPM_RSEM.rds')
TPM = TPM[rowSums(TPM >= 5)> 37/4,]
eDat = log2(TPM+1)


library(WGCNA)
# choose a soft threshold power
powers = c(seq(from = 1, to=31, by=2));
sft=pickSoftThreshold(t(eDat),dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")

# Plot the results
pdf('results/WGCNA_SoftThresh.pdf', height=5, width=10)
sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Generating adjacency and TOM similarity matrices based on the selected softpower
softPower = 29

# building WGCNA modules
net = blockwiseModules(t(eDat), power = 29, #0.8SFT
                       TOMType = "signed", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       saveTOMFileBase = "DataTOM",
                       verbose = 3)
#saveRDS(net, 'results/net_gene_logTPM.rds')
net = readRDS('results/net_gene_logTPM.rds')
dynamicMods = net$colors
dynamicColors = paste0('module ',dynamicMods)
dynamicColors = recode(dynamicColors, `module 0`='module 00',
                       `module 1`='module 01',`module 2`='module 02',
                       `module 3`='module 03',`module 4`='module 04',
                       `module 5`='module 05',`module 6`='module 06',
                       `module 7`='module 07',`module 8`='module 08',
                       `module 9`='module 09')
table(dynamicColors)

# log2TPM
#module 00 module 01 module 02 module 03 module 04 module 05 module 06 module 07
#    12380      1831      1069      1042       840       321       203       197
#module 08 module 09 module 10 module 11 module 12 module 13
#      186       127        41        35        23        19

# extract modules
gene.names=rownames(eDat)
module_colors= unique(dynamicColors)
module=list()
for (color in module_colors){
  module[[color]]=gene.names[which(dynamicColors==color)]
}
mod=list()
for (color in module_colors){
  mod[[color]]=as.data.frame(cbind(color,gene.names[which(dynamicColors==color)]))
}

# Quantify module similarity by eigengene correlation. Eigengenes: Module representatives
MEList = moduleEigengenes(t(eDat), colors = dynamicColors)
MEs = MEList$eigengenes
MEs = orderMEs(MEs)

png('results/eigen_dis.png',height=5000, width=5000, res=600)
plotEigengeneNetworks(MEs, "", marDendro = c(0,6,1,4), marHeatmap = c(5,6,1,2))
dev.off()

# cell-specific markers
MEs = MEs %>% dplyr::select(-`MEmodule 00`)
design <- model.matrix(~0+MType+Method+SMOKING_STATUS, data=pDat)
corfit <- duplicateCorrelation(t(MEs),design,block=pDat$SampleName)
corfit$consensus #intra-subject correlation
fit <- lmFit(t(MEs),design,block=pDat$SampleName,correlation=corfit$consensus)
cm <- makeContrasts(
		MTypeDN-(MTypeDP+MTypeM1+MTypeM2)/3,
		MTypeDP-(MTypeDN+MTypeM1+MTypeM2)/3,
		MTypeM1-(MTypeDN+MTypeDP+MTypeM2)/3,
		MTypeM2-(MTypeDN+MTypeDP+MTypeM1)/3,
		levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

topt_DN = topTable(fit2, coef=1, num=Inf) %>% mutate(probe_id=rownames(.))
fwrite(topt_DN, 'results/WGCNA_DN_logTPM.csv')

topt_DP = topTable(fit2, coef=2, num=Inf) %>% mutate(probe_id=rownames(.))
fwrite(topt_DP, 'results/WGCNA_DP_logTPM.csv')

topt_M1 = topTable(fit2, coef=3, num=Inf) %>% mutate(probe_id=rownames(.))
fwrite(topt_M1, 'results/WGCNA_M1_logTPM.csv')

topt_M2 = topTable(fit2, coef=4, num=Inf) %>% mutate(probe_id=rownames(.))
fwrite(topt_M2, 'results/WGCNA_M2_logTPM.csv')


# combine the results into a table of estimates and a table of p-values
library(purrr)
sigmodule = rbind(topt_DN %>% filter(adj.P.Val<0.1),topt_DP %>% filter(adj.P.Val<0.1),
				  topt_M1 %>% filter(adj.P.Val<0.1),topt_M2 %>% filter(adj.P.Val<0.1)) %>% dplyr::select(probe_id) %>%
				  unique() %>% unlist()
Est = purrr::reduce(list(
				topt_DN %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, logFC),
				topt_DP %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, logFC),
				topt_M1 %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, logFC),
				topt_M2 %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, logFC)),
				function(x,y) merge(x,y, by='probe_id'))
rownames(Est) = Est$probe_id
colnames(Est) = c('Probe','DN','DP','M1','M2')
Est$Probe=NULL

FDR = purrr::reduce(list(
				topt_DN %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, adj.P.Val),
				topt_DP %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, adj.P.Val),
				topt_M1 %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, adj.P.Val),
				topt_M2 %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, adj.P.Val)),
				function(x,y) merge(x,y, by='probe_id'))
rownames(FDR) = FDR$probe_id
colnames(FDR) = c('Probe','DN','DP','M1','M2')
FDR$Probe=NULL

# make a modules vs phenotypes heatmap
color_dis = table(dynamicColors) %>% as.data.frame() %>% filter(dynamicColors %in% gsub('ME','',sigmodule))
color_dis$dynamicColors=as.character(color_dis$dynamicColors)
color_dis$dynamicColors == gsub('ME','',rownames(Est))
lab_dis = paste0(color_dis$dynamicColors,' (',color_dis$Freq, ')')

png('results/module_marker_logTPM.png',height=6000, width=4500, res=600)
# Will display the estimates and their p-values
textMatrix = paste(signif(as.matrix(Est), 2), "\n(",
                   signif(as.matrix(FDR), 3), ")", sep = "");
dim(textMatrix) = dim(Est)
par(mar = c(6, 10, 3,2));
labeledHeatmap(Matrix = as.matrix(Est),
               xLabels = c('DN','DP','M1','M2'),
               yLabels = lab_dis,
               ySymbols = lab_dis,
               colorLabels = FALSE,
               cex.lab.x = 1.2,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.8,
               zlim = c(-0.15,0.15),
               main = paste(""))
dev.off()

# module enrichment
mod_genes = rbindlist(mod)
mod_genes = mod_genes %>% mutate(Gene = anno$gene_name[match(V2,anno$gene_id)])
ref = mod_genes$Gene %>% .[!is.na(.) & .!=''] %>% unique()
mod_enrich_res = lapply(gsub('ME','',sigmodule[-1]), function(x) mod_genes %>% filter(color %in% x) %>% dplyr::select(Gene) %>% unlist() %>% .[!is.na(.) & .!=''] %>% unique() %>% enrich.analysis(.,ref))
names(mod_enrich_res) = gsub('ME','',sigmodule[-1])
fwrite(mod_genes,'results/module_genes_logTPM.csv')

mod_enrich_GOBP=list()
for(i in 1:length(mod_enrich_res)){
	mod_enrich_GOBP[[i]] = mod_enrich_res[[i]]$GOBP %>% simplify() %>% .@result %>% dplyr::filter(p.adjust<0.1) %>%
					mutate(module = names(mod_enrich_res)[i])
}
fwrite(mod_enrich_GOBP %>% rbindlist(), 'results/WGCNA_module_enrichment_GOBP_logTPM_nRundt.csv')

mod_enrich_KEGG=list()
for(i in 1:length(mod_enrich_res)){
	mod_enrich_KEGG[[i]] = mod_enrich_res[[i]]$KEGG@result %>% dplyr::filter(p.adjust<0.1) %>%
					mutate(module = names(mod_enrich_res)[i])
}
fwrite(mod_enrich_KEGG %>% rbindlist(), 'results/WGCNA_module_enrichment_KEGG_logTPM.csv')

# plot top GOBP pathways in each module
dat=fread('results/WGCNA_module_enrichment_GOBP_logTPM_nRundt.csv')
plotdat = as.data.table(dat %>% arrange(module))[,head(.SD %>% filter(Count>2) %>% arrange(p.adjust),3),by=module]
plotdat = plotdat %>% mutate(number=1:nrow(plotdat))
setnames(plotdat,'p.adjust','FDR')
ggplot(plotdat, aes(x=module, y=reorder(Description,number), size=Count, fill=-log10(FDR)))+
  geom_point(shape=21)+
  scale_size_area(name = 'Number of Overlapping Genes', max_size = 10)+
  scale_fill_gradient(low = "#e8d3ff", high = "#240049", space = "Lab",
                      na.value = "grey50", guide = "colourbar",  name = expression(paste('- Log'[10],' FDR'))) +
  xlab('Modules')+
  ylab('GO Biological Processes')+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(axis.title=element_text(face="bold",size=11),axis.text=element_text(size=11,face="bold"))
ggsave('results/enrichment_GOBP_modulemarker.png',height=7, width = 15)

# WGCNA MM
# names (colors) of the modules
modNames = substring(names(MEs), 3)
# calculation
geneModuleMembership = as.data.frame(cor(t(eDat), MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),37));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
write.csv(geneModuleMembership,'MMmembership_num_logTPM.csv')
write.csv(MMPvalue,'MMmembershipP_num_logTPM.csv')

all(rownames(geneModuleMembership)==rownames(MMPvalue))
MM=cbind(geneModuleMembership %>% mutate(Gene=rownames(.)) %>% dplyr::select(Gene,everything()),
         MMPvalue)
MMlis=list()
for(i in 2:14){
  MMlis[[i-1]]=MM[,c(1,i,i+13)][order(MM[,i],decreasing = T),]
}
# top 20 MM
do.call(cbind,lapply(MMlis, function(x) head(x,20) %>% mutate(Symbol = anno$gene_name[match(Gene, anno$gene_id)], `EmptyCol`=NA))) %>% fwrite(.,'top20modgene_num_logTPM.csv')

## plot module genes by pathway
design <- model.matrix(~0+MType+Method+SMOKING_STATUS.x, data=pDat)
corfit <- duplicateCorrelation(eDat,design,block=pDat$SampleName)
corfit$consensus #intra-subject correlation
fit <- lmFit(eDat,design,block=pDat$SampleName,correlation=corfit$consensus)
fit2 <- eBayes(fit)
dat=coefficients(fit2)[,1:4]

plotGivenGenes = function(dt, genes, height, width, prefix){
	genes2Plot = anno[match(genes, anno$gene_name),]
	exprs = dt[genes2Plot$gene_id,]
	rownames(exprs)=genes2Plot$gene_name

	input_final = t(scale(t(exprs)))
	rownames(input_final)=genes2Plot$gene_name

	trim_threshold = 2  # trim exp
	input_final[input_final>trim_threshold] = trim_threshold
	input_final[input_final< (-trim_threshold)] = -trim_threshold

	aheatmap(input_final,Rowv =NA, Colv = NA, labCol = c('DN','DP','M1','M2'),  border_color='white',
		 filename = paste0("results/logTPM/",prefix,"_heatmap.png"),
         width = width, height = height)
}

mtg=c('NDUFAB1/UQCRC1/NDUFB4/NDUFB2/ATP5F1D/NDUFB7/NDUFC1/NDUFS8/COX6A1/NDUFS7/COX7A2L/PARK7/ATP5PB/NDUFB3/NDUFA8/ATP5F1E/NDUFA1/COX6B1/COX7C/UQCR11/NDUFA5/COX4I1/COX7B/NDUFA2/COX5B/UQCC2/NDUFA9/NDUFB10/NDUFS6/NIPSNAP2/NDUFB9/SURF1/NDUFC2/ATP5F1A/ATP5MC3/ATP5PF/UQCRB/NDUFS2/ATP5MC1/NDUFS4/UQCRQ/COX6C/NDUFB6/STOML2/ATP5F1C/NDUFB8/ATP5MG/NDUFV1/ATP5PD/SDHAF2/NDUFS5/ATP5ME/UQCRFS1/NDUFA3/CYCS/UQCRH/NDUFA11/COX8A/NDUFV2/COX5A/CYC1/NDUFB1/UQCR10/NDUFA12/NDUFA13/NDUFA4/NDUFS3/ATP5MF/ATP5PO/NDUFA7')
mtg = c('NDUFAB1/UQCRC1/NDUFB4/NDUFB2/ATP5F1D/NDUFB7/NDUFC1/NDUFS8/COX6A1/NDUFS7/COX7A2L/PARK7/ATP5PB/NDUFB3/NDUFA8/ATP5F1E/NDUFA1/COX6B1/COX7C/UQCR11/NDUFA5/COX4I1/COX7B/NDUFA2/COX5B/UQCC2/NDUFA9/NDUFB10/NDUFS6/NIPSNAP2/NDUFB9/SURF1/NDUFC2/ATP5F1A/ATP5MC3/')
mtg=c('ATP5PF/UQCRB/NDUFS2/ATP5MC1/NDUFS4/UQCRQ/COX6C/NDUFB6/STOML2/ATP5F1C/NDUFB8/ATP5MG/NDUFV1/ATP5PD/SDHAF2/NDUFS5/ATP5ME/UQCRFS1/NDUFA3/CYCS/UQCRH/NDUFA11/COX8A/NDUFV2/COX5A/CYC1/NDUFB1/UQCR10/NDUFA12/NDUFA13/NDUFA4/NDUFS3/ATP5MF/ATP5PO/NDUFA7')
mtg = strsplit(mtg,split='/') %>% unlist()
plotGivenGenes(dat, mtg, 6.5, 3, 'M02_OxidativePhosphorylation_2')
cyt = c('PIK3CB/CYFIP2/ELMO2/ABL1/MAPK1/PIK3R2/WASL/WIPF1/PIK3CA/NCKAP1L/ACTR2/PIK3R1/DOCK1/ELMO1/WASF2/PRKCD/SYK/WIPF2/SRC/PLCG2')
cyt = strsplit(cyt,split='/') %>% unlist()
plotGivenGenes(dat, cyt, 6.5, 3, 'M01_phagocytosis_test')
myl = c('PSMB1/MVP/SERPINB1/ERP44/CYBA/GDI2/ACAA1/AP1M1/GSTP1/DYNLL1/CST3/VAPA/PGRMC1/PSMD7/COTL1/PSMA2/PSMD3/LAMTOR3/COMMD9/NIT2/RAP1A/LAMTOR2/CD58/VAMP8/NPC2/YPEL5/CYSTM1/SERPINB6/RAP1B/GMFG/LGALS3/OSTF1/CD63/RAC1/DBNL/PSMB7/SRP14/PFKL/PSMA5/ILF2/DEGS1/ARL8A/STK11IP/DYNLT1/PTGES2/COMMD3/LAMTOR1/DSN1/NDUFC2/FCER1G/ARPC5/S100A11/PSMD6/GYG1/APEH/BRI3/PSMC3/B2M/RAB4B/AGPAT2/KCMF1/DPP7/IST1/ANXA2/PSMD13/TUBB4B/HMGB1/PPIA/APRT/CSNK2B/MIF/NME2')
myl = strsplit(myl,split='/') %>% unlist()
plotGivenGenes(dat, myl, 6.5, 3, 'M02_myeloidCellActivation')
ret = c('NDUFAB1/UQCRC1/NDUFB4/NDUFB2/NDUFB7/ETFB/NDUFC1/NDUFS8/COX6A1/NDUFS7/COX7A2L/PARK7/SDHB/NDUFB3/NDUFA8/NDUFA1/COX6B1/COX7C/UQCR11/NDUFA5/COX4I1/COX7B/NDUFA2/COX5B/NDUFA9/NDUFB10/NDUFS6/NDUFB9/NDUFC2/UQCRB/NDUFS2/NDUFS4/UQCRQ/COX6C/NDUFB6/NDUFB8/NDUFV1/SDHAF2/NDUFS5/UQCRFS1/NDUFA3/CYCS/UQCRH/NDUFA11/COX8A/NDUFV2/COX5A/CYC1/NDUFB1/UQCR10/NDUFA12/NDUFA13/NDUFA4/NDUFS3/NDUFA7')
ret = strsplit(ret,split='/') %>% unlist()
plotGivenGenes(dat, ret, 6.5, 3, 'M02_RespElectronTransport')
fao = c('PDK4/CROT/ALDH3A2/ACACB/ECHDC1/MLYCD/DECR1/ECH1/SCP2/ACADM/ACADS/POR/HACL1/PPARG/PRKAA1/HADHB/HADH/AKT1/TYSND1/ETFDH/IRS2')
fao = strsplit(fao,split='/') %>% unlist()
plotGivenGenes(dat, fao, 6.5, 3, 'M03_FattyAcidOxidation')
cc = c('AKAP8L/ANLN/TACC3/ANAPC4/ASPM/AURKA/TPX2/ANAPC5/BIRC5/KIF4A/MYBL2/CDC25B/AKAP8/NCAPG/CUL9/SMC4/WRAP73/CDC20/CENPF/KIF14/ZWINT/PKMYT1/TUBGCP6/SGO1/TOP2A/CCNB1/NUSAP1/KIF11/CENPE/KIF2C/FANCD2/CDCA5/MKI67/BUB1B/MAD2L1/CDT1/BUB1/MUS81/SH2B1/AURKB/L3MBTL1/KIF18B/SIRT7/EME2/PRC1/SPIRE2/KIFC1')
cc = strsplit(cc,split='/') %>% unlist()
plotGivenGenes(dat, cc, 6.5, 3, 'M04_NuclearDivision')
tm = c('CCT4/TCP1/CCT7/CCT6A/CCT2')
tm = strsplit(tm,split='/') %>% unlist()
plotGivenGenes(dat, tm, 3, 3, 'M06_ProteinLocalizationToTelemere')
sg=c('FOS/HSPA1A/JUN/FOSB/JUNB/EGR1/HSPA1B/UBC/ZFP36/HSPB1/HSP90AA1/MT2/DNAJB1/BTG2/NR4A1/CEBPD/HSPA8/MT1/LER2/DNAJA1/SOCS3/ATF3/JUND/CEBPD/LD3/PPP1R15A/HSPE1/CXCL1/DUSP1/HSP90AB1/NFKBIA/HSPH1')
sg = strsplit(sg,split='/') %>% unlist() %>% .[. %in% anno$gene_name]
plotGivenGenes(dat, sg, 6.5, 3, 'stress_gene')







############################## additional analysis ##################
# DE on CD163- vs CD163+
pDat = pDat %>% mutate(CD163=ifelse(MType %in% c('M2','DP'),'Positive','Negative'))
design <- model.matrix(~CD163+Method+SMOKING_STATUS.x, data=pDat)
corfit <- duplicateCorrelation(eDat,design,block=pDat$SampleName)
corfit$consensus #intra-subject correlation
fit <- lmFit(eDat,design,block=pDat$SampleName,correlation=corfit$consensus)
fit2 <- eBayes(fit)

topt = topTable(fit2, coef=2, num=Inf)
topt_annot = topt %>% mutate(Gene_ID=rownames(.)) %>% mutate(Gene=anno$gene_name[match(Gene_ID,anno$gene_id)])
fwrite(topt_annot, 'results/CD163_neg_vs_pos.csv')

DEgenes = topt_annot %>% filter(adj.P.Val<0.1)  #128

pcutoff=topt_annot %>% filter(adj.P.Val<0.1) %>% dplyr::select(P.Value) %>% unlist() %>% max()
labcutoff=topt_annot %>% filter(logFC>0) %>% head(20) %>% dplyr::select(P.Value) %>% unlist() %>% max()
p=volcano_plot_exprs(topt_annot,
					contrast.fld = 'logFC', pval.fld = 'P.Value',
                    pthresh = pcutoff, contrast.thresh = 0,
                    xlabel = expression(paste(log[2], ' Fold Change')),
                    ylabel = expression(paste('-',log[10], ' P')),
                    label.fld = 'Gene',
                    line.label='FDR=0.1',
                    label.thresh = labcutoff, label.thresh.fld = 'P.Value')+
	theme(legend.position='none')
ggsave('results/volcano_CD163_neg_vs_pos2.png', height =5, width=5)

# WGCNA on CD163- vs CD163+
net = readRDS('results/net_gene_logTPM.rds')
dynamicMods = net$colors
dynamicColors = paste0('module ',dynamicMods)
dynamicColors = recode(dynamicColors, `module 0`='module 00',
                       `module 1`='module 01',`module 2`='module 02',
                       `module 3`='module 03',`module 4`='module 04',
                       `module 5`='module 05',`module 6`='module 06',
                       `module 7`='module 07',`module 8`='module 08',
                       `module 9`='module 09')
table(dynamicColors)

# log2TPM
#module 00 module 01 module 02 module 03 module 04 module 05 module 06 module 07
#    12380      1831      1069      1042       840       321       203       197
#module 08 module 09 module 10 module 11 module 12 module 13
#      186       127        41        35        23        19

# extract modules
gene.names=rownames(eDat)
module_colors= unique(dynamicColors)
module=list()
for (color in module_colors){
  module[[color]]=gene.names[which(dynamicColors==color)]
}
mod=list()
for (color in module_colors){
  mod[[color]]=as.data.frame(cbind(color,gene.names[which(dynamicColors==color)]))
}

# Quantify module similarity by eigengene correlation. Eigengenes: Module representatives
MEList = moduleEigengenes(t(eDat), colors = dynamicColors)
MEs = MEList$eigengenes
MEs = orderMEs(MEs)

# cell-specific markers
MEs = MEs %>% dplyr::select(-`MEmodule 00`)
design <- model.matrix(~CD163+Method+SMOKING_STATUS.x, data=pDat)
corfit <- duplicateCorrelation(t(MEs),design,block=pDat$SampleName)
corfit$consensus #intra-subject correlation
fit <- lmFit(t(MEs),design,block=pDat$SampleName,correlation=corfit$consensus)
fit2 <- eBayes(fit)

topt = topTable(fit2, coef=2, num=Inf) %>% mutate(probe_id=rownames(.)) %>% arrange(probe_id)
fwrite(topt, 'results/WGCNA_CD163_logTPM.csv')

# combine the results into a table of estimates and a table of p-values
library(purrr)
sigmodule = rbind(topt %>% filter(adj.P.Val<0.1)) %>% dplyr::select(probe_id) %>%
				  unique() %>% unlist()
Est = purrr::reduce(list(
				topt %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, logFC)),
				function(x,y) merge(x,y, by='probe_id'))
rownames(Est) = Est$probe_id
colnames(Est) = c('Probe','CD163+')
Est$Probe=NULL

FDR = purrr::reduce(list(
				topt %>% filter(probe_id %in% sigmodule) %>% dplyr::select(probe_id, adj.P.Val)),
				function(x,y) merge(x,y, by='probe_id'))
rownames(FDR) = FDR$probe_id
colnames(FDR) = c('Probe','CD163+')
FDR$Probe=NULL

# make a modules vs phenotypes heatmap
color_dis = table(dynamicColors) %>% as.data.frame() %>% filter(dynamicColors %in% gsub('ME','',sigmodule))
color_dis$dynamicColors=as.character(color_dis$dynamicColors)
color_dis$dynamicColors == gsub('ME','',rownames(Est))
lab_dis = paste0(color_dis$dynamicColors,' (',color_dis$Freq, ')')

png('results/module_marker_logTPM_CD163.png',height=3000, width=3000, res=600)
# Will display the estimates and their p-values
textMatrix = paste(signif(as.matrix(Est), 2), "\n(",
                   signif(as.matrix(FDR), 3), ")", sep = "");
dim(textMatrix) = dim(Est)
par(mar = c(6, 10, 3,2));
labeledHeatmap(Matrix = as.matrix(Est),
               xLabels = c('CD163+'),
               yLabels = lab_dis,
               ySymbols = lab_dis,
               colorLabels = FALSE,
               cex.lab.x = 1.2,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.8,
               zlim = c(-0.15,0.15),
               main = paste(""))
dev.off()


# expression clustering based on IL6, IL1B, ARG1, KLF4
genes=c('IL6','IL1B','ARG1','KLF4','CD40','CD163')
anno[match(genes,anno$gene_name)]
pDat = pDat %>% mutate(#IL6=as.numeric(eDat['ENSG00000136244.12_5',]),
					   IL1B=as.numeric(eDat['ENSG00000125538.11_5',]),
					   #ARG1=as.numeric(eDat['ENSG00000118520.14_2',]),
					   KLF4=as.numeric(eDat['ENSG00000136826.15_5',]),
					   CD40=as.numeric(eDat['ENSG00000101017.13_6',]),
					   CD163=as.numeric(eDat['ENSG00000177575.12_3',]))
ggplot(pDat,aes(x=IL1B,y=KLF4,color=MType))+
  geom_point()+
  scale_color_manual(name = 'Macrophage Subtype', values = c('blue','orange','red','purple')) +
  theme_bw()+
  theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'),
          legend.position = 'top',
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))+
  stat_ellipse(level=0.68,type='norm')
ggsave('results/macrophage_by_IL1B_KLF4.png', height =5, width=5)

aheatmap(eDat[c('ENSG00000125538.11_5','ENSG00000136826.15_5'),],labCol = pDat$MType,
  border_color='white',filename = "results/addMarkers_heatmap.png",
  width = 7, height = 5)

wd <- "~/Desktop/ohta_add_20191223/20191223"
setwd(wd)
Sys.setenv('R_MAX_VSIZE'=80000000000)


library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)


#データの読み込み
file<-"GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"

#１行目と２行目が変なので３行目から読み込む
ex<-data.table::fread(file, stringsAsFactors=FALSE, sep="\t", data.table = FALSE,skip=2)
rownames(ex)<-ex[,1]
e<-ex[,-c(1,16293)]

#１行目を読み込む
L1<-readLines(con = file, n = 1)
L1<-unlist(strsplit(L1,"\t"))

#２行目を読み込む
L2<-readLines(con = file, n = 2)
L2<-unlist(strsplit(L2,"\t"))

#発現値データに列名、行名をつける
colnames(e)<-L1[-1]

#サンプル情報の読み込み
df<-read.csv("GSE120575_patient_ID_single_cells.txt",skip=19,sep="\t",as.is=T)
df<-df[,!is.na(df[1,])]

#データの整理
df<-df[1:16291,]
rownames(df)<-df[,2]
common<-intersect(colnames(e),rownames(df))
df<-df[common,]
e<-e[,common]
df[,1]<-sub(" ","",df[,1])
colnames(e)<-df[,1]
rownames(df)<-df[,1]

#Seuratオブジェクトの作成
EX <- CreateSeuratObject(counts = e, project = "GSE120575", min.cells = 0)
rm(e)


#ミトコンドリア遺伝子の含有率を調べる
EX[["percent.mt"]] <- PercentageFeatureSet(EX, pattern = "^MT-")

#遺伝子数、カウント数、ミトコンドリア遺伝子含有率を可視化
#png("VlnPlot.png")
VlnPlot(EX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
#dev.off()

#正規化ここは飛ばしていいはず
#Supplementary_files_format_and_content: tab-delimited text file containing log2(TPM+1) values (transcript-per-million reads) for 55,737 transcripts (rows) and 16,291 cells (columns)
#EX <- NormalizeData(EX, normalization.method = "LogNormalize", scale.factor = 10000)

#分散を調べる
EX <- FindVariableFeatures(EX, selection.method = "vst", nfeatures = 2000)

#スケーリングする
all.genes <- rownames(EX)
EX <- ScaleData(EX, features = all.genes)

#主成分スコアの計算
EX <- RunPCA(EX, features = VariableFeatures(object = EX))

#クラスターを探す
EX <- FindNeighbors(EX, dims = 1:20)
EX <- FindClusters(EX, resolution = 0.8)

#UMAPの計算
EX <- RunUMAP(EX, dims = 1:20)
#png("UMAP.png")
DimPlot(EX, reduction = "umap")
dev.off()

##複数遺伝子のラベル付け
gene1 <- "TCF7" #commandArgs(trailingOnly = TRUE)[1]
gene2 <- "PDCD1" #commandArgs(trailingOnly = TRUE)[2]
#method <- "median" #commandArgs(trailingOnly = TRUE)[3] #implemented methods are median/mean/zero

#カウントデータ
e2 <- as.matrix(EX@assays$RNA@counts)

gene1_median_rm0 <- median(e2[gene1,][(e2[gene1,])!=0])
gene2_median_rm0 <- median(e2[gene2,][(e2[gene2,])!=0])

newgroup <- factor(case_when(
e2[gene1,] >= gene1_median_rm0 & e2[gene2,] >= gene2_median_rm0 ~ "High_High",
e2[gene1,] > 0 & e2[gene1,] < gene1_median_rm0 & e2[gene2,] >= gene2_median_rm0 ~ "Low_High",
e2[gene1,] == 0 & e2[gene2,] >= gene2_median_rm0 ~ "Minus_High",
e2[gene1,] >= gene1_median_rm0 & e2[gene2,] > 0 & e2[gene2,] < gene2_median_rm0 ~ "High_Low",
e2[gene1,] > 0 & e2[gene1,] < gene1_median_rm0 & e2[gene2,] > 0 & e2[gene2,] < gene2_median_rm0 ~ "Low_Low",
e2[gene1,] == 0 & e2[gene2,] > 0 & e2[gene2,] < gene2_median_rm0 ~ "Minus_Low",
e2[gene1,] >= gene1_median_rm0 & e2[gene2,] == 0 ~ "High_Minus",
e2[gene1,] > 0 & e2[gene1,] < gene1_median_rm0 & e2[gene2,] == 0 ~ "Low_Minus",
e2[gene1,] == 0 & e2[gene2,] ==0 ~ "Minus_Minus"))
                      
EX[["newgroup"]] <- newgroup

#複数遺伝子のラベル付け
#png("UMAP_split.png")
DimPlot(EX, reduction = "umap", group.by = "newgroup")
dev.off()

#DEG
un_g <- unique(newgroup)

for (i in (1:length(un_g))){
tem_group <- ifelse(newgroup==un_g[i], as.character(un_g[i]),"others")

EX[["tem_group"]] <- tem_group

genes.deg <- FindMarkers(EX, ident.1 =as.character(un_g[i]),ident.2 = "others", group.by = "tem_group")
genes.deg$avg_logFC <- 2^genes.deg$avg_logFC
colnames(genes.deg) [2] <- "FC"

colnames(genes.deg)[3] <- paste0(as.character(un_g[i]),"_ratio")
colnames(genes.deg)[4] <-"Other_ratio"

upDEG <- genes.deg[which(genes.deg$FC > 1 & genes.deg$p_val_adj < 0.05),]
downDEG <- genes.deg[which(genes.deg$FC < 1 & genes.deg$p_val_adj < 0.05),]

out_f1 <- paste0(as.character(un_g[i]), "DEG.csv")
out_f2 <- paste0(as.character(un_g[i]), "upDEG.csv")
out_f3 <- paste0(as.character(un_g[i]), "downDEG.csv")

setwd(wd)
setwd("./2gene_result/")
dir.create(as.character(un_g[i]))
setwd(as.character(un_g[i]))
write.csv(genes.deg, out_f1)
write.csv(upDEG, out_f2)
write.csv(downDEG, out_f3)
setwd(wd)
}



##########
#ここまで

genes.deg <- FindMarkers(EX, ident.1 = "DoublePositive", ident.2 = "DoubleNegative",  group.by = "newgroup")
genes.deg$avg_logFC <- 2^genes.deg$avg_logFC
colnames(genes.deg) [2] <- "FC"
write.csv(deg, "DP_DN.csv")


#1 MM_vs_HM

#2 MM_vs_LM

#3 MM_vs_ML

#4 MM_vs_LL

#5 MM_vs_HL

#6 MM_vs_MH

#7 MM_vs_LH

#8 MM_vs_HH



#9 HM_vs_LM

#10 HM_vs_ML

#11 HM_vs_LL

#12 HM_vs_HL

#13 HM_vs_MH

#14 HM_vs_LH

#15 HM_vs_HH


#LM_vs_ML

#LM_vs_LL

#




#responderラベル指定
responder.group <- factor(df[,6])
EX[["responder.group"]] <- responder.group
#responderのDEG探し
genes_responder.deg <- FindMarkers(EX, ident.1 = "Responder", ident.2 = "Non-responder",  group.by = "responder.group")
genes_responder.deg$avg_logFC <- 2^genes_responder.deg$avg_logFC
colnames(genes_responder.deg) [2] <- "FC"
write.csv(genes_responder.deg, "responder_DEG.csv")


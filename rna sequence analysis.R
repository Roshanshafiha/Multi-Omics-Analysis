#load the required library

library("limma")
library("dplyr")
library("ggplot2")
library("ggrepel")
library( "EnhancedVolcano")
library( "AnnotationDbi" )
library("org.Hs.eg.db")
library("biomaRt")
library("tibble")

#control versus 1x concentration 

homosapien_annotation<-read.csv("ortho_dma_to_hsa_mappings.csv",header = F)

count_data<- read.csv(file = "rna_vst_counts.csv")

rownames(count_data)<-count_data$X

count_data<-count_data[,-c(1)]

metadata<-read.csv("sample_sheet.csv")

#remove the 10x from the count and metadata 

metadata_1x<-metadata[metadata$REF != "10x", ] 

count_data_1x<-subset(count_data, select=metadata_1x$SampleID)

#Get the differentially expressed genes

#create the design matrix

design1 <- model.matrix(~0+metadata_1x$REF,ref="Control")

#rename the columns

colnames(design1) <- c("normal","one")

fit_data1 <- lmFit(count_data_1x, design1)


contrasts1 <- makeContrasts(normal - one, levels=design1)

fit2_data1 <- contrasts.fit(fit_data1, contrasts1)

fit2_data1<- eBayes(fit2_data1)

full_results1 <- topTable(fit2_data1, number=Inf)

full_results1 <- cbind(symbols = rownames(full_results1), full_results1)

rownames(full_results1) <- 1:nrow(full_results1)


#get the rownames of human and match it with daphnia 

full_results1<-cbind(homosapien_protein= homosapien_annotation$V2[ match(full_results1$symbols,homosapien_annotation$V1) ],
      full_results1)


#annotate the plot according to the ensemble (Homosapien genes)

p_cutoff <- 0.05
fc_cutoff <- 1
topN <-50

full_results1 %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN,homosapien_protein,"")) %>% 
  ggplot(aes(x = logFC, y = B , col=Significant,label=Label)) + geom_point() +
          geom_text_repel(col="black")

#extract the significant genes with proper annotations 

sign_genes1<-full_results1 %>% filter(full_results1$adj.P.Val < 0.05 & 
                                        abs(full_results1$logFC) > 1)

sign_genes1<- sign_genes1[complete.cases(sign_genes1), ]

protein_ID_1x<-sign_genes1$homosapien_protein


# define biomart object

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# get the genes for the given protein IDs

results_1x <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","ensembl_peptide_id"),
                 filters = "ensembl_peptide_id", values = protein_ID_1x,
                 mart = mart)

#match the protein names with the gene names.

sign_genes1<-cbind(homosapien_gene= results_1x$ensembl_gene_id[ match(sign_genes1$homosapien_protein,results_1x$ensembl_peptide_id) ],
                     sign_genes1)


anno <- AnnotationDbi::select(org.Hs.eg.db, sign_genes1$homosapien_gene,columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), keytype="ENSEMBL")

sign_genes1$gene_ID <- anno$SYMBOL

sign_genes1<-sign_genes1 %>% add_column(gene_ID = anno$SYMBOL, .before = 1)

write.csv(sign_genes1,"1x samples human annotations.csv")




#heatmap generation

#extract the genes greater then 0.05

selected  <- sign_genes1$symbols

expression_data_heatmap <- count_data_1x[ rownames(count_data_1x) %in% selected, ]

heatmap(as.matrix(expression_data_heatmap))





#control versus 10x concentration

#remove the 1x from the count and metadata 

metadata_10x<-metadata[metadata$REF != "1x", ] 

count_data_10x<-subset(count_data, select=metadata_10x$SampleID)

#Get the differentially expressed genes

#create the design matrix

design2 <- model.matrix(~0+metadata_10x$REF,ref="Control")

#rename the columns

colnames(design2) <- c("normal","ten")

fit_data2 <- lmFit(count_data_10x, design2)


contrasts2 <- makeContrasts(normal - ten, levels=design2)

fit2_data2 <- contrasts.fit(fit_data2, contrasts2)

fit2_data2<- eBayes(fit2_data2)

full_results2 <- topTable(fit2_data2, number=Inf)

full_results2 <- cbind(symbols = rownames(full_results2), full_results2)

rownames(full_results2) <- 1:nrow(full_results2)


#get the rownames of human and match it with daphnia 

full_results2<-cbind(homosapien_protein= homosapien_annotation$V2[ match(full_results2$symbols, homosapien_annotation$V1) ],
                     full_results2)


#annotate the plot according to the ensemble (Homosapien genes)

p_cutoff <- 0.05
fc_cutoff <- 1
topN <-50

full_results2 %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN,homosapien_protein,"")) %>% 
  ggplot(aes(x = logFC, y = B , col=Significant,label=Label)) + geom_point() +
  geom_text_repel(col="black")



sign_genes2<-full_results2 %>% filter(full_results2$adj.P.Val < 0.05 & 
                                        abs(full_results2$logFC) > 1)

sign_genes2<- sign_genes2[complete.cases(sign_genes2), ]

protein_ID_2x<-sign_genes2$homosapien_protein


# define biomart object

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# get the genes for the given protein IDs

results_2x <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id","ensembl_peptide_id"),
                    filters = "ensembl_peptide_id", values = protein_ID_2x,
                    mart = mart)

#match the protein names with the gene names.

sign_genes2<-cbind(homosapien_gene= results_2x$ensembl_gene_id[ match(sign_genes2$homosapien_protein,results_2x$ensembl_peptide_id) ],
                   sign_genes2)


anno <- AnnotationDbi::select(org.Hs.eg.db, sign_genes2$homosapien_gene,columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), keytype="ENSEMBL")

sign_genes2$gene_ID <- anno$SYMBOL

sign_genes2<-sign_genes2 %>% add_column(gene_ID = anno$SYMBOL, .before = 1)

write.csv(sign_genes2,"10x samples human annotations.csv")


#heatmap

selected_10x <- sign_genes2$symbols

expression_data_heatmap_10x <- count_data_10x[ rownames(count_data_10x) %in% selected_10x, ]

heatmap(as.matrix(expression_data_heatmap_10x))






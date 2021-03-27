# Multi-Omics-Analysis
The advanced machine learning approaches and their integrative analysis for a complex real-world study

# Dataset provided :

River organic extracts:

A table of 38 chemicals on the row and 12 sites on the column and the concentration values at environmental level (1x) in the cell (10x level could be calculated based on these data).

Transcriptome profiling:

A total of 150 samples were successfully sequenced, including 6 control samples, 72 1x-treatment samples, and 72 10x-treatment samples. All the reads were mapped to transcriptome reference database via Salmon, further normalised by DESeq2 package. The data table was with 18,903 genes on the row, 150 sample on the column and the normalised read count number in the cell.

Metabolome profiling:

positive mode, a data table with 1,285 peaks on the row, 149 samples on the column with their processed peak intensity values in the cell.
negative mode, a data table with 2,331 peaks on the row, 145 samples on the column with processed peak intensity values.

Omics feature annotation:

The Daphnia magna genes are mapped to the other model species for more convenient annotation. The ortholog mapping is done using the latest version of OrthoDB to the Daphnia pulex , Drosophila melanogaster and Homo sapiens. 

Research Goals:

Answer the following questions.

1.How do the inter-connected biological features (transcriptome and/or metabolome) change between the different locations (D01-D12), concentration levels (Control, 1x, 10x), and their interactions?
2.How does the change relate to the distribution of detected individual organic chemical compound? Which chemical may cause the most significant adverse outcome on the Daphnia?
3.What biological/environmental insights can you obtain from your data analysis findings? 

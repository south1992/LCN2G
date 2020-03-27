# LC-N2G(LC-Vis): A Local Consistency Approach for Nutrigenomics Data Analysis

    This shiny app was created for Visualization of nutrigenomics data analysis. The shiny app composed of 3 parts: 
Data preprocessing, Gene clustering and Gene-Nutrition visualization. 
    * Data preprocessing: This part works for preprocess the gene and nutrition data. In the left penal,
one can choose the filter threshold for both mean and standard deviant of gene expression as well as the absolute Spearman correlation between gene and nutrition.
The normalization method can also be choosed as cpm or log2 cpm. Output density plot of mean and variance of gene expression and absolute Spearman correlation were shown in right panel/
Under the right panel, the number of gene after filtering is shown.
    * Gene Clustering: In order to get representitive gene for next step visualization, we provided a gene cluster method. If we have background know

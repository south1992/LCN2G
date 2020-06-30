# LC-N2G: A Local Consistency Approach for Nutrigenomics Data Analysis


    It's a shiny app(LC-N2G) to explore the relationship of nutrition and its corresponding gene expression data. The default dataset come from mouse nutrition study(GSE85998). The overall workflow of LC-N2G works as follows:
	
	![LC-Vis](/img/fig1.png)
	
	<div align=center>
	Figure 1 Overall workflow of LC-N2G. Gn×p2 and Nn×p1 represent input of matrix of gene and nutrition respectively. First step we calculate LC-Stat of combinations with a gene of interest to find combination of nutrients with small LC-Stat. Then a LC-Test is performed to evaluate the relationship between combination of nutrients with gene. Finally the NGF is performed for selected combination and genes.
	</div>
	
## Packages Requirements

- shiny ≥ 1.4.0.2
- shinythemes ≥ 1.1.2
- shinyjs ≥ 1.1
- DT ≥ 0.13
- directlabels ≥ 2020.1.31
- mclust ≥ 5.4.5
- GA ≥ 3.2
- WGCNA ≥ 1.69
- tidyverse ≥ 1.3.0
- plotly ≥ 4.9.2
- visNetwork ≥ 2.0.9
- limma ≥ 3.42.2
- dynamicTreeCut ≥ 1.63-1
- psych ≥ 1.9.12.31
- fields ≥ 10.3
- plotly ≥ 4.9.2
- reshape2 ≥ 1.4.3

 
## Usage	

	This shiny app is composed in 3 parts, Data Preprocessing, Gene Clustering and Visualization. We introduce each part as follows.

### Data Preprocessing

	This part works for preprocess the gene and nutrition data. 
 
	![LC-Vis](/img/git_fig2.png)
	
	We have implemented different kind of preprocess method in this part. One can upload its own dataset(in csv format and each row for a sample) or use the default dataset which upload nothing and click Analysis(GSE85998). Different threshold for filtering is provided and the resulting plot will appear in the right panel after click analysis.
	In the default dataset, we alreadly normalised using Cyclophilin as the endogenous control, so just use "None" as normalization method.

### Gene Clustering

    In order to get representitive gene for next step visualization, we provided a gene cluster method. If some special gene expression are interested, one can go to next step.

	![LC-Vis](/img/git_fig3.png)	
	
	The clustering method we used here are hierachical clustering and the choices of distance matrix varies. The first panel provided a clustering method for gene-gene relationship and the distance matrix can based on gene expression or gene-gene correlation. Also the gene-nutrition relationship can use to clustering either based on covariance or correlation matrix. A hybrid parameter was used to combine this two matrix. If we only use one, we can set this parameter either be 0 or 1. The cuttree method also can be choosed.
	The output panel will visualize the dendrogram for clustering and the corresponding result can be find in tabel tab-panel. A network of nutrition, clusters and its corresponding eigen-gene is presented in summary tab-panel.
	
### Gene Nutrition Visualization

	Finally, the nutrition geometry framework(NGF) is done in this part.
	
	![LC-Vis](/img/git_fig4.png)
	
	The first parameter here are used for choose gene or cluster of gene. If a gene is used, the result will be NGF of the gene or if a cluster choosed, NGFs of mean, variance, eigen-gene and a presentitive gene will visualize.
	For LC-Test of the relationship between selected gene/cluster, the candidate nutrition variables can be choosed in the second panel. Click analysis in this panel, the combination with smallest LC-Statistic will show in the output panel.
	The final parameter are used to visualize the NGF with selected xaxis and yaxis. Output panel will show the resulting NGF.

---
title: "Introduction to Single Cell RNAseq Part 5"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: March 24 2021, 5pm

# Part 5: 

## Load libraries

```r
library(Seurat)
library(ggplot2)
library(dplyr)
```

## Load the Seurat object

```r
load(file="pca_sample_corrected.RData")
experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 36601 features across 4000 samples within 1 assay 
 Active assay: RNA (36601 features, 2000 variable features)
  1 dimensional reduction calculated: pca
</div>

## So how many features should we use? Use too few and your leaving out interesting variation that may define cell types, use too many and you add in noise? maybe?

Lets choose the first 50, based on our prior part.


```r
use.pcs = 1:50
```

## Identifying clusters

Seurat implements an graph-based clustering approach. Distances between the cells are calculated based on previously identified PCs. 

The default method for identifying k-nearest neighbors has been changed in V4 to [annoy](https://github.com/spotify/annoy) ("Approximate Nearest Neighbors Oh Yeah!). This is an approximate nearest-neighbor approach that is widely used for high-dimensional analysis in many fields, including single-cell analysis. Extensive community benchmarking has shown that annoy substantially improves the speed and memory requirements of neighbor discovery, with negligible impact to downstream results. 



Seurat prior approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNAseq data. Briefly, Seurat identified clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors (KNN) and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. You can switch back to using the previous default setting using nn.method="rann".


The FindClusters function implements the neighbor based clustering procedure, and contains a resolution parameter that sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters. I tend to like to perform a series of resolutions, investigate and choose.


```r
?FindNeighbors
```


```r
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction="pca", dims = use.pcs)

experiment.aggregate <- FindClusters(
    object = experiment.aggregate,
    resolution = seq(0.25,4,0.5),
    verbose = FALSE
)
```


Seurat add the clustering information to the metadata beginning with RNA_snn_res. followed by the resolution


```r
head(experiment.aggregate[[]])
```

<div class='r_output'>                        orig.ident nCount_RNA nFeature_RNA batchid percent.mito
 CCGTTCAAGGTGACCA-PBMC2      PBMC2       1533         1004  Batch1     4.044357
 ATTGGTGTCGGTTCGG-PBMC2      PBMC2       1225          834  Batch1     7.673469
 GTAGGCCTCTATCCTA-PBMC2      PBMC2       1109          806  Batch1     4.418395
 CATTATCCAGGCAGTA-PBMC2      PBMC2       1949         1226  Batch1     7.439713
 CTGGTCTTCCAAATGC-PBMC2      PBMC2       3074         1303  Batch1     2.374756
 AGCTTGAGTTGTGGAG-PBMC2      PBMC2       1028          748  Batch1     3.599222
                             S.Score   G2M.Score Phase old.ident
 CCGTTCAAGGTGACCA-PBMC2 -0.009757111  0.05114394   G2M     PBMC2
 ATTGGTGTCGGTTCGG-PBMC2 -0.006619531  0.02301644   G2M     PBMC2
 GTAGGCCTCTATCCTA-PBMC2  0.019239276 -0.02686039     S     PBMC2
 CATTATCCAGGCAGTA-PBMC2 -0.077384425  0.04775564   G2M     PBMC2
 CTGGTCTTCCAAATGC-PBMC2  0.032593746 -0.01541937     S     PBMC2
 AGCTTGAGTTGTGGAG-PBMC2 -0.051954726 -0.01295370    G1     PBMC2
                        RNA_snn_res.0.25 RNA_snn_res.0.75 RNA_snn_res.1.25
 CCGTTCAAGGTGACCA-PBMC2                3                3                5
 ATTGGTGTCGGTTCGG-PBMC2                3                3                5
 GTAGGCCTCTATCCTA-PBMC2                2               10                9
 CATTATCCAGGCAGTA-PBMC2                3                3                5
 CTGGTCTTCCAAATGC-PBMC2                2                6                6
 AGCTTGAGTTGTGGAG-PBMC2                9               12               13
                        RNA_snn_res.1.75 RNA_snn_res.2.25 RNA_snn_res.2.75
 CCGTTCAAGGTGACCA-PBMC2                5                5                3
 ATTGGTGTCGGTTCGG-PBMC2                5                5                3
 GTAGGCCTCTATCCTA-PBMC2                9                8                8
 CATTATCCAGGCAGTA-PBMC2                5                5                3
 CTGGTCTTCCAAATGC-PBMC2                6                6                6
 AGCTTGAGTTGTGGAG-PBMC2               13               15               18
                        RNA_snn_res.3.25 RNA_snn_res.3.75 seurat_clusters
 CCGTTCAAGGTGACCA-PBMC2                4                6               6
 ATTGGTGTCGGTTCGG-PBMC2                4                6               6
 GTAGGCCTCTATCCTA-PBMC2                8                8               8
 CATTATCCAGGCAGTA-PBMC2                4                6               6
 CTGGTCTTCCAAATGC-PBMC2                6                5               5
 AGCTTGAGTTGTGGAG-PBMC2               18               18              18
</div>

Lets first investigate how many clusters each resolution produces and set it to the smallest resolutions of 0.5 (fewest clusters).


```r
sapply(grep("res",colnames(experiment.aggregate@meta.data),value = TRUE),
       function(x) length(unique(experiment.aggregate@meta.data[,x])))
```

<div class='r_output'> RNA_snn_res.0.25 RNA_snn_res.0.75 RNA_snn_res.1.25 RNA_snn_res.1.75 
               10               15               18               23 
 RNA_snn_res.2.25 RNA_snn_res.2.75 RNA_snn_res.3.25 RNA_snn_res.3.75 
               25               27               27               30
</div>
### Plot TSNE coloring for each resolution

tSNE dimensionality reduction plots are then used to visualize clustering results. As input to the tSNE, you should use the same PCs as input to the clustering analysis.


```r
experiment.aggregate <- RunTSNE(
  object = experiment.aggregate,
  reduction.use = "pca",
  dims.use = use.pcs,
  do.fast = TRUE)
```



```r
DimPlot(object = experiment.aggregate, group.by=grep("res",colnames(experiment.aggregate@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/tsne_all_1-1.png)<!-- -->


```r
DimPlot(object = experiment.aggregate, group.by=grep("res",colnames(experiment.aggregate@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/tsne_all_2-1.png)<!-- -->

1. Try exploring different PCS, so first 5, 15, 25, we used 50, what about 100? How does the clustering change?

Once complete go back to 1:50

### Choosing a resolution

Lets set the default identity to a resolution of 0.25 and produce a table of cluster to sample assignments.

```r
Idents(experiment.aggregate) <- "RNA_snn_res.0.25"
table(Idents(experiment.aggregate),experiment.aggregate$orig.ident)
```

<div class='r_output'>    
     PBMC2 PBMC3 T021PBMC T022PBMC
   0    37    18      438      527
   1    58   709       12        9
   2   429    12        0        3
   3   307    42       32       13
   4     9     5      295       67
   5     7     7       86      232
   6    77    14       91       29
   7    22   165        6       17
   8     0     4       36      100
   9    54    24        4        3
</div>
Plot TSNE coloring by the slot 'ident' (default).

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "tsne", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/plot_tsne-1.png)<!-- -->


### uMAP dimensionality reduction plot.


```r
experiment.aggregate <- RunUMAP(
  object = experiment.aggregate,
  dims = use.pcs)
```

Plot uMap coloring by the slot 'ident' (default).

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, reduction = "umap", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/plot_umap-1.png)<!-- -->

Catagorical data can be plotted using the DimPlot function.


TSNE plot by cell cycle

```r
DimPlot(object = experiment.aggregate, pt.size=0.5, group.by = "Phase", reduction = "umap" )
```

![](scRNA_Workshop-PART5_files/figure-html/plot_cellcycle-1.png)<!-- -->

1. Try creating a table, of cluster x cell cycle

### Can use feature plot to plot our read valued metadata, like nUMI, Feature count, and percent Mito
FeaturePlot can be used to color cells with a 'feature', non categorical data, like number of UMIs

```r
FeaturePlot(experiment.aggregate, features = c('nCount_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/plot_rna-1.png)<!-- -->
and number of genes present

```r
FeaturePlot(experiment.aggregate, features = c('nFeature_RNA'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/plot_feature-1.png)<!-- -->

percent mitochondrial

```r
FeaturePlot(experiment.aggregate, features = c('percent.mito'), pt.size=0.5)
```

![](scRNA_Workshop-PART5_files/figure-html/plot_mito-1.png)<!-- -->

## Building a phylogenetic tree relating the 'average' cell from each group in default 'Ident' (currently "RNA_snn_res.0.25"). Tree is estimated based on a distance matrix constructed in either gene expression space or PCA space.


```r
experiment.aggregate <- BuildClusterTree(
  experiment.aggregate, dims = use.pcs)

PlotClusterTree(experiment.aggregate)
```

![](scRNA_Workshop-PART5_files/figure-html/create_tree-1.png)<!-- -->
1. Create new trees of other data

Once complete go back to Res 0.25


```r
DimPlot(object = experiment.aggregate, pt.size=0.5, label = TRUE, reduction = "umap")
```

![](scRNA_Workshop-PART5_files/figure-html/umap_plot2-1.png)<!-- -->

### Merging clusters

Merge Clustering results, so lets say clusters 3 and 7 are actually the same cell type (overlaps in the tsne, not as apparent in the umap) and we don't wish to separate them out as distinct clusters. Same with 4 and 5.


```r
experiment.merged = experiment.aggregate
Idents(experiment.merged) <- "RNA_snn_res.0.25"

experiment.merged <- RenameIdents(
  object = experiment.merged,
  '3' = '7', '4' = '5'
)

table(Idents(experiment.merged))
```

<div class='r_output'> 
    7    5    0    1    2    6    8    9 
  604  708 1020  788  444  211  140   85
</div>
```r
DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "umap")
```

![](scRNA_Workshop-PART5_files/figure-html/merging_cluster-1.png)<!-- -->

```r
VlnPlot(object = experiment.merged, features = "percent.mito", pt.size = 0.05)
```

![](scRNA_Workshop-PART5_files/figure-html/merging_cluster-2.png)<!-- -->

### Reording the clusters

In order to reorder the clusters for plotting purposes take a look at the levels of the Ident, which indicates the ordering, then relevel as desired.


```r
experiment.examples <- experiment.merged
levels(experiment.examples@active.ident)
```

<div class='r_output'> [1] "7" "5" "0" "1" "2" "6" "8" "9"
</div>
```r
experiment.examples@active.ident <- relevel(experiment.examples@active.ident, "6")
levels(experiment.examples@active.ident)
```

<div class='r_output'> [1] "6" "7" "5" "0" "1" "2" "8" "9"
</div>
```r
# now cluster 6 is the "first" factor

DimPlot(object = experiment.merged, pt.size=0.5, label = T, reduction = "umap")
```

![](scRNA_Workshop-PART5_files/figure-html/merging_cluster2-1.png)<!-- -->

```r
VlnPlot(object = experiment.merged, features = "percent.mito", pt.size = 0.05)
```

![](scRNA_Workshop-PART5_files/figure-html/merging_cluster2-2.png)<!-- -->



```r
# relevel all the factors to the order I want
Idents(experiment.examples) <- factor(experiment.examples@active.ident, levels=c("6","2","0","7","5","9","1","8"))
levels(experiment.examples@active.ident)
```

<div class='r_output'> [1] "6" "2" "0" "7" "5" "9" "1" "8"
</div>
```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "umap")
```

![](scRNA_Workshop-PART5_files/figure-html/merging_cluster3-1.png)<!-- -->


### Re-assign clustering result (subclustering only cluster 0) to clustering for resolution 3.75  (@ reslution 0.25) [adding a R prefix]

```r
newIdent = as.character(Idents(experiment.examples))
newIdent[newIdent == '0'] = paste0("R",as.character(experiment.examples$RNA_snn_res.3.75[newIdent == '0']))

Idents(experiment.examples) <- as.factor(newIdent)
table(Idents(experiment.examples))
```

<div class='r_output'> 
   1   2   5   6   7   8   9 R12 R15 R18 R19  R2 R24  R3  R8  R9 
 788 444 708 211 604 140  85 132 101   2  44 270  46 256   5 164
</div>

```r
DimPlot(object = experiment.examples, pt.size=0.5, label = T, reduction = "umap")
```

![](scRNA_Workshop-PART5_files/figure-html/subclusters_plot-1.png)<!-- -->

Plot UMAP  coloring by the slot 'orig.ident' (sample names) with alpha colors turned on. A pretty picture

```r
DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "umap" )
```

![](scRNA_Workshop-PART5_files/figure-html/pretty_pre-1.png)<!-- -->



```r
## Pretty umap using alpha
alpha.use <- 2/5
p <- DimPlot(object = experiment.aggregate, group.by="orig.ident", pt.size=0.5, reduction = "umap")
p$layers[[1]]$mapping$alpha <- alpha.use
p + scale_alpha_continuous(range = alpha.use, guide = F)
```

![](scRNA_Workshop-PART5_files/figure-html/pretty_post-1.png)<!-- -->

Removing cells assigned to clusters from a plot, So here plot all clusters but cluster 6 (contaminant?)

```r
# create a new tmp object with those removed
experiment.aggregate.tmp <- experiment.aggregate[,-which(Idents(experiment.aggregate) %in% c("6"))]

dim(experiment.aggregate)
```

<div class='r_output'> [1] 36601  4000
</div>
```r
dim(experiment.aggregate.tmp)
```

<div class='r_output'> [1] 36601  3789
</div>


```r
DimPlot(object = experiment.aggregate.tmp, group.by="orig.ident", pt.size=0.5, reduction = "umap", label = T)
```

![](scRNA_Workshop-PART5_files/figure-html/removing_cells_plot-1.png)<!-- -->

## Identifying Marker Genes

Seurat can help you find markers that define clusters via differential expression.

`FindMarkers` identifies markers for a cluster relative to all other clusters.

`FindAllMarkers` does so for all clusters

`FindAllMarkersNode` defines all markers that split a Node from the cluster tree


```r
?FindMarkers
```


```r
markers = FindMarkers(experiment.merged, ident.1=c(6))

head(markers)
```

<div class='r_output'>       p_val avg_log2FC pct.1 pct.2 p_val_adj
 FCRL1     0   2.042489 0.498 0.004         0
 FCRLA     0   1.476939 0.389 0.001         0
 AFF3      0   2.329102 0.682 0.018         0
 BANK1     0   2.416954 0.711 0.009         0
 PAX5      0   1.492776 0.455 0.003         0
 MS4A1     0   3.771599 0.910 0.014         0
</div>
```r
dim(markers)
```

<div class='r_output'> [1] 2029    5
</div>
```r
table(markers$avg_log2FC > 0)
```

<div class='r_output'> 
 FALSE  TRUE 
  1360   669
</div>
```r
table(markers$p_val_adj < 0.05)
```

<div class='r_output'> 
 FALSE  TRUE 
  1016  1013
</div>

pct.1 and pct.2 are the proportion of cells with expression above 0 in ident.1 and ident.2 respectively. p_val is the raw p_value associated with the differntial expression test with adjusted value in p_val_adj. avg_logFC is the average log fold change difference between the two groups.

avg_diff (lines 130, 193 and) appears to be the difference in log(x = mean(x = exp(x = x) - 1) + 1) between groups.  It doesn’t seem like this should work out to be the signed ratio of pct.1 to pct.2 so I must be missing something.  It doesn’t seem to be related at all to how the p-values are calculated so maybe it doesn’t matter so much, and the sign is probably going to be pretty robust to how expression is measured.

Can use a violin plot to visualize the expression pattern of some markers

```r
VlnPlot(object = experiment.merged, features = rownames(markers)[1:2], pt.size = 0.05)
```

![](scRNA_Workshop-PART5_files/figure-html/vln-1.png)<!-- -->

Or a feature plot

```r
FeaturePlot(
    experiment.merged,
    head(rownames(markers), n=2),
    cols = c("lightgrey", "blue"),
    ncol = 2
)
```

![](scRNA_Workshop-PART5_files/figure-html/gene_feature-1.png)<!-- -->

FindAllMarkers can be used to automate the process across all genes.


```r
markers_all <- FindAllMarkers(
    object = experiment.merged,
    only.pos = TRUE,
    min.pct = 0.25,
    thresh.use = 0.25
)
dim(markers_all)
```

<div class='r_output'> [1] 4398    7
</div>
```r
head(markers_all)
```

<div class='r_output'>               p_val avg_log2FC pct.1 pct.2     p_val_adj cluster  gene
 CTSW  1.556588e-233   2.001720 0.957 0.344 5.697267e-229       7  CTSW
 GZMA  1.132526e-188   1.715444 0.896 0.302 4.145157e-184       7  GZMA
 PRF1  6.009217e-186   1.812446 0.829 0.235 2.199433e-181       7  PRF1
 NKG7  4.692792e-174   1.532782 0.978 0.401 1.717609e-169       7  NKG7
 KLRD1 8.593663e-166   1.565966 0.724 0.186 3.145367e-161       7 KLRD1
 CST7  2.633715e-154   1.379715 0.868 0.290 9.639660e-150       7  CST7
</div>
```r
table(table(markers_all$gene))
```

<div class='r_output'> 
    1    2    3    4    5 
 1441 1063  254   16    1
</div>
```r
markers_all_single <- markers_all[markers_all$gene %in% names(table(markers_all$gene))[table(markers_all$gene) == 1],]

dim(markers_all_single)
```

<div class='r_output'> [1] 1441    7
</div>
```r
table(table(markers_all_single$gene))
```

<div class='r_output'> 
    1 
 1441
</div>
```r
table(markers_all_single$cluster)
```

<div class='r_output'> 
   7   5   0   1   2   6   8   9 
 145 187 165 303  44 158 189 250
</div>
```r
head(markers_all_single)
```

<div class='r_output'>               p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
 GZMK   2.404916e-82  1.5066984 0.298 0.055 8.802235e-78       7   GZMK
 CLSTN3 7.169670e-49  0.8599396 0.255 0.069 2.624171e-44       7 CLSTN3
 MYBL1  3.189092e-48  0.9405144 0.276 0.082 1.167240e-43       7  MYBL1
 ATG2A  4.199373e-43  0.8775682 0.318 0.115 1.537013e-38       7  ATG2A
 PIK3R1 1.509989e-40  0.9409462 0.675 0.464 5.526710e-36       7 PIK3R1
 CD8B   2.363168e-34  0.7106495 0.270 0.092 8.649433e-30       7   CD8B
</div>
Plot a heatmap of genes by cluster for the top 10 marker genes per cluster

```r
top10 <- markers_all_single %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(
    object = experiment.merged,
    features = top10$gene
)
```

![](scRNA_Workshop-PART5_files/figure-html/markers_head-1.png)<!-- -->


```r
# Get expression of genes for cells in and out of each cluster
getGeneClusterMeans <- function(gene, cluster){
  x <- GetAssayData(experiment.merged)[gene,]
  m <- tapply(x, ifelse(Idents(experiment.merged) == cluster, 1, 0), mean)
  mean.in.cluster <- m[2]
  mean.out.of.cluster <- m[1]
  return(list(mean.in.cluster = mean.in.cluster, mean.out.of.cluster = mean.out.of.cluster))
}

## for sake of time only using first six (head)
means <- mapply(getGeneClusterMeans, head(markers_all[,"gene"]), head(markers_all[,"cluster"]))
means <- matrix(unlist(means), ncol = 2, byrow = T)

colnames(means) <- c("mean.in.cluster", "mean.out.of.cluster")
rownames(means) <- head(markers_all[,"gene"])
markers_all2 <- cbind(head(markers_all), means)
head(markers_all2)
```

<div class='r_output'>               p_val avg_log2FC pct.1 pct.2     p_val_adj cluster  gene
 CTSW  1.556588e-233   2.001720 0.957 0.344 5.697267e-229       7  CTSW
 GZMA  1.132526e-188   1.715444 0.896 0.302 4.145157e-184       7  GZMA
 PRF1  6.009217e-186   1.812446 0.829 0.235 2.199433e-181       7  PRF1
 NKG7  4.692792e-174   1.532782 0.978 0.401 1.717609e-169       7  NKG7
 KLRD1 8.593663e-166   1.565966 0.724 0.186 3.145367e-161       7 KLRD1
 CST7  2.633715e-154   1.379715 0.868 0.290 9.639660e-150       7  CST7
       mean.in.cluster mean.out.of.cluster
 CTSW         2.774018           0.7832086
 GZMA         2.277050           0.6690290
 PRF1         2.249473           0.5942009
 NKG7         3.743689           1.2028493
 KLRD1        1.470524           0.3646510
 CST7         2.307280           0.7399939
</div>
## Finishing up clusters.

At this point in time you should use the tree, markers, domain knowledge, and goals to finalize your clusters. This may mean adjusting PCA to use, mergers clusters together, choosing a new resolutions, etc. When finished you can further name it cluster by something more informative. Ex.

```r
experiment.clusters <- experiment.aggregate
experiment.clusters <- RenameIdents(
  object = experiment.clusters,
  '0' = 'cell_type_A',
  '1' = 'cell_type_B',
  '2' = 'cell_type_C'
)
# and so on

DimPlot(object = experiment.clusters, pt.size=0.5, label = T, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/finish_cluster-1.png)<!-- -->

Right now our results ONLY exist in the Ident data object, lets save it to our metadata table so we don't accidentally loose it.

```r
experiment.merged$finalcluster <- Idents(experiment.merged)
```

## Subsetting samples and plotting

If you want to look at the representation of just one sample, or sets of samples

```r
experiment.sample2 <- subset(experiment.merged, orig.ident == "T021PBMC")

DimPlot(object = experiment.sample2, group.by = "RNA_snn_res.0.25", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/subset-1.png)<!-- -->

```r
experiment.batch1 <- subset(experiment.merged, batchid == "Batch1")

DimPlot(object = experiment.batch1, group.by = "RNA_snn_res.0.25", pt.size=0.5, label = TRUE, reduction = "tsne")
```

![](scRNA_Workshop-PART5_files/figure-html/subset-2.png)<!-- -->

### Adding in a new metadata column representing samples within clusters. So differential expression of PBMC2 vs PBMC3 within cluster 7


```r
experiment.merged$samplecluster = paste(experiment.merged$orig.ident,experiment.merged$finalcluster,sep = '-')

# set the identity to the new variable
Idents(experiment.merged) <- "samplecluster"

markers.comp <- FindMarkers(experiment.merged, ident.1 = "PBMC2-7", ident.2= "PBMC3-7")

head(markers.comp)
```

<div class='r_output'>               p_val avg_log2FC pct.1 pct.2    p_val_adj
 RPS4Y1 2.287891e-74  -2.824205 0.000 0.749 8.373911e-70
 MT-CO3 2.417890e-62   1.550807 1.000 0.981 8.849721e-58
 AREG   9.795578e-59  -3.270285 0.088 0.768 3.585280e-54
 IFITM1 7.799421e-56  -1.229551 0.954 0.995 2.854666e-51
 LY6E   3.730028e-49  -1.678001 0.611 0.952 1.365227e-44
 IFITM2 1.457371e-43  -1.234761 0.824 0.990 5.334125e-39
</div>
```r
experiment.subset <- subset(experiment.merged, samplecluster %in%  c( "PBMC2-7", "PBMC3-7" ))
DoHeatmap(experiment.subset, features = head(rownames(markers.comp),20))
```

![](scRNA_Workshop-PART5_files/figure-html/unnamed-chunk-1-1.png)<!-- -->



```r
Idents(experiment.merged) <- "finalcluster"
```

And last lets save all the objects in our session.

```r
save(list=ls(), file="clusters_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-March-Single-Cell-RNA-Seq-Analysis/master/data_analysis/scRNA_Workshop-PART6.Rmd", "scRNA_Workshop-PART6.Rmd")
```

## Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.0.3 (2020-10-10)
 Platform: x86_64-apple-darwin17.0 (64-bit)
 Running under: macOS Big Sur 10.16
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
 LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base     
 
 other attached packages:
 [1] dplyr_1.0.5        ggplot2_3.3.3      SeuratObject_4.0.0 Seurat_4.0.1      
 
 loaded via a namespace (and not attached):
   [1] Rtsne_0.15            colorspace_2.0-0      deldir_0.2-10        
   [4] ellipsis_0.3.1        ggridges_0.5.3        spatstat.data_2.1-0  
   [7] leiden_0.3.7          listenv_0.8.0         farver_2.1.0         
  [10] ggrepel_0.9.1         RSpectra_0.16-0       fansi_0.4.2          
  [13] codetools_0.2-18      splines_4.0.3         knitr_1.31           
  [16] polyclip_1.10-0       jsonlite_1.7.2        ica_1.0-2            
  [19] cluster_2.1.1         png_0.1-7             uwot_0.1.10          
  [22] shiny_1.6.0           sctransform_0.3.2     spatstat.sparse_2.0-0
  [25] compiler_4.0.3        httr_1.4.2            assertthat_0.2.1     
  [28] Matrix_1.3-2          fastmap_1.1.0         lazyeval_0.2.2       
  [31] limma_3.44.3          later_1.1.0.1         htmltools_0.5.1.1    
  [34] tools_4.0.3           igraph_1.2.6          gtable_0.3.0         
  [37] glue_1.4.2            RANN_2.6.1            reshape2_1.4.4       
  [40] Rcpp_1.0.6            scattermore_0.7       jquerylib_0.1.3      
  [43] vctrs_0.3.6           ape_5.4-1             nlme_3.1-152         
  [46] lmtest_0.9-38         xfun_0.22             stringr_1.4.0        
  [49] globals_0.14.0        mime_0.10             miniUI_0.1.1.1       
  [52] lifecycle_1.0.0       irlba_2.3.3           goftest_1.2-2        
  [55] future_1.21.0         MASS_7.3-53.1         zoo_1.8-9            
  [58] scales_1.1.1          spatstat.core_1.65-5  promises_1.2.0.1     
  [61] spatstat.utils_2.1-0  parallel_4.0.3        RColorBrewer_1.1-2   
  [64] yaml_2.2.1            reticulate_1.18       pbapply_1.4-3        
  [67] gridExtra_2.3         sass_0.3.1            rpart_4.1-15         
  [70] stringi_1.5.3         highr_0.8             rlang_0.4.10         
  [73] pkgconfig_2.0.3       matrixStats_0.58.0    evaluate_0.14        
  [76] lattice_0.20-41       ROCR_1.0-11           purrr_0.3.4          
  [79] tensor_1.5            patchwork_1.1.1       htmlwidgets_1.5.3    
  [82] labeling_0.4.2        cowplot_1.1.1         tidyselect_1.1.0     
  [85] parallelly_1.24.0     RcppAnnoy_0.0.18      plyr_1.8.6           
  [88] magrittr_2.0.1        R6_2.5.0              generics_0.1.0       
  [91] DBI_1.1.1             pillar_1.5.1          withr_2.4.1          
  [94] mgcv_1.8-34           fitdistrplus_1.1-3    survival_3.2-10      
  [97] abind_1.4-5           tibble_3.1.0          future.apply_1.7.0   
 [100] crayon_1.4.1          KernSmooth_2.23-18    utf8_1.2.1           
 [103] spatstat.geom_2.0-1   plotly_4.9.3          rmarkdown_2.7        
 [106] grid_4.0.3            data.table_1.14.0     digest_0.6.27        
 [109] xtable_1.8-4          tidyr_1.1.3           httpuv_1.5.5         
 [112] munsell_0.5.0         viridisLite_0.3.0     bslib_0.2.4
</div>
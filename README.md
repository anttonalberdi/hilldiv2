# HillDiv 2

hilldiv2 is an R package that provides a set of functions to assist analysis of diversity for diet reconstruction, microbial community profiling or more general ecosystem characterisation analyses based on Hill numbers, using OTU/ASV/MAG tables, and associated phylogenetic trees and trait/attribute tables as inputs. The package includes functions for measurement of neutral, phylogenetic and functional Hill numbers, partitioning of neutral, phylogenetic and functional Hill numbers, and dissimilarity measurements derived from them.

The statistical framework developed around Hill numbers encompasses many of the most broadly employed diversity (e.g. richness, Shannon index, Simpson index), phylogenetic diversity (e.g. Faith’s PD, Allen’s H, Rao’s quadratic entropy) and dissimilarity (e.g. Sørensen index, Unifrac distances) metrics. This enables the most common analyses of diversity to be performed while grounded in a single statistical framework. For details about the use of Hill numbers in molecularly characterised biological systems, you can read the following article:

Alberdi A, Gilbert MTP. (2019). A guide to the application of Hill numbers to DNA‐based diversity analyses. Molecular Ecology Resources. 19(4): 804-817. https://doi.org/10.1111/1755-0998.13014

For relevant literature visit the [References](#references)

## Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Data input](#data-input)
  - [Hill numbers diversity metrics using **hilldiv()**](#hill-numbers-diversity-metrics)
  - [Hill numbers diversity partitioning using **hillpart()**](#hill-numbers-diversity-partitioning)
  - [Hill numbers similarity metrics using **hillsim()**](#hill-numbers-similarity-metrics)
  - [Hill numbers dissimilarity metrics using **hilldiss()**](#hill-numbers-dissimilarity-metrics)
  - [Hill numbers pairwise dissimilarities using **hillpart()**](#hill-numbers-pairwise-dissimilarities)
  - [Hill numbers redundancy metrics using **hillred()**](#hill-numbers-redundancy-metrics)
- [References](#references)

## Installation

hilldiv2 can be installed from this Github repository using devtools.

```r
install.packages("devtools")
library(devtools)
install_github("anttonalberdi/hilldiv2")
library(hilldiv2)
library(tidyverse)
```

## Usage

The following examples use the lizard-associated metagenome-assembled genome (MAG) dataset included in the package.

### Data input

#### Counts table

The counts table needs to be in a tabular format, with OTUs/ASVs/MAGs in rows, and samples in columns. Note that hilldiv2 only does total sum squares (TSS) normalisation to ensure the data become proportional. In some contexts, the data might require some further filtering or normalisation for reliable estimation of diversity.

|      | Sample1 | Sample2 | Sample3 | Sample4 |
| ---  |---|---|---|---|
| MAG1 | 342   | 24522 | 0     | 1232  |
| MAG2 | 0     | 0     | 32424 | 732   |
| MAG3 | 82242 | 22342 |232    | 52232 |
| MAG4 |223    | 1214  |12113  | 0     |


#### Phylogenetic tree
The phylogenetic tree needs to be in "phylo" format. This requires a newick- or nexus-format tree to be loaded to R using the library **ape**.

#### Trait table

The trait table needs to be in a tabular format, with OTUs/ASVs/MAGs in rows, and traits in columns. The traits can be continuous, binary or proportional, but the appropriate metric for generating the distance matrix will need to be chosen.

|      | Trait1 | Trait2 | Trait3 | Trait4 |
| ---  |---|---|---|---|
| MAG1 | 1   | 0.2 | 0.7  | 1    |
| MAG2 | 0   | 0   | 0.2  | 0.98 |
| MAG3 | 0.2 | 0   | 0.13 | 0.80 |
| MAG4 | 0.9 | 1   | 0.8  | 0.7  |

#### Sample metadata table

The sample metadata table needs to be in a tabular format, with samples in rows, and relevant metadata fields in columns.

|      | Population | Sex | Weight | Length |
| ---  |---|---|---|---|
| Sample1 | PopulationA | Female | 7.54  | 12.4 |
| Sample2 | PopulationA | Male   | 6.24  | 11.2 |
| Sample3 | PopulationB | Female | 8.43  | 14.2 |
| Sample4 | PopulationB | Male   | 4.35  | 10.1 |

### Load data and prepare functional data

```r
#Load data
counts <- lizards_counts
tree <- lizards_tree
traits <- lizards_traits
metadata <- lizards_metadata

#Convert traits into distance matrix
dist <- traits2dist(traits, method="gower")
```

### Hill numbers diversity metrics

The versatile function ***hilldiv()*** serves to compute Hill numbers of different types – neutral, phylogenetic, and functional – based on the input provided. When only count data is supplied, the function calculates neutral Hill numbers, which are also known as taxonomic Hill numbers. When both count data and a phylogenetic tree are provided, the function computes phylogenetic Hill numbers, in the form of "effective number of lineages" or qD(T) according to Chao et al. (2016). If a distance matrix is included as input, the function calculates functional Hill numbers. In the absence of a specific q-value, the function defaults to calculating Hill numbers of order q=0, q=1, and q=2 to measure diversity.

```r
hilldiv(data=counts)
hilldiv(data=counts,tree=tree)
hilldiv(data=counts,dist=dist)
```

If the order(s) of diversity is/are defined, then Hill numbers are computed for the desired q-values.

```r
hilldiv(data=counts,q=0)
hilldiv(data=counts,q=c(0,0.5,1),tree=tree)
hilldiv(data=counts,q=2,dist=dist)
```

For functional Hill numbers, it is possible to define a desired Tau value as well. The default Tau value is the maximum distance between any two OTUs/ASVs/MAGs.

```r
hilldiv(data=counts,q=2,dist=dist,tau=0.3)
hilldiv(data=counts,q=2,dist=dist,tau=max(dist))
```

#### Hill numbers diversity plotting

Hill numbers can be then displayed visually using ggplot2.

```r
hilldiv(data=counts,q=1,tree=tree) %>%
    as.data.frame() %>%
    pivot_longer(everything(),names_to = "sample", values_to = "value") %>%
    left_join(metadata, by = join_by(sample == sample)) %>%
    ggplot(., aes(x=population, y=value, color=population, fill=population)) +
        geom_boxplot() +
        geom_jitter() +
        scale_color_manual(values = c("#E3D97B","#46edc8","#374d7c")) +
        scale_fill_manual(values = c("#E3D97B50","#46edc850","#374d7c50")) +
        labs(y = "Phylogenetic diversity") +
        theme_classic() +
        theme(legend.position="null",
        axis.title.x = element_blank())
```

![Visualisation of phylogenetic diversities of q=1](/images/alpha.png)

### Hill numbers diversity partitioning

The function ***hillpart()*** performs partitioning of any type of Hill numbers-based diversity (neutral, phylogenetic and functional) at any order(s) of diversity (q>0) into alpha, gamma and beta components.

```r
hillpart(data=counts)
hillpart(data=counts,tree=tree)
hillpart(data=counts,dist=dist)

hillpart(data=counts,q=0)
hillpart(data=counts,q=c(0,0.5,1),tree=tree)
hillpart(data=counts,q=2,dist=dist)
```

### Hill numbers similarity metrics

The function ***hillsim()*** computes overall similarity metrics based on the Hill numbers beta diversity following Chiu et al. (2014). By default the function outputs the following four similarity metrics, although users can limit the output to a single or some of those metrics. Calculations can be done for a single q-value or multiple q-values. Note these metrics are the complements (1-X) of the dissimilarity metrics computed by ***hilldiss()***.

- **S** ***(Jaccard-type turnover-complement):*** it quantifies the complement of the normalised species turnover rate in a sample relative to the total pool of samples (i.e., gamma diversity).
- **V** ***(Sørensen-type turnover-complement):*** it quantifies the complement of the normalised species turnover rate in a sample relative to one sample (i.e., alpha diversity).
- **U** ***(Jaccard-type overlap):*** it quantifies the proportion of shared species in the total pool of samples. Therefore, this metric quantifies similarity from the perspective of the pool of samples.
- **C** ***(Sørensen-type overlap):*** it quantifies the effective average proportion of shared OTUs/ASVs/MAGs in samples. Therefore, this metric quantifies similarity from the perspective of a single sample.

```r
hillsim(data=counts)
hillsim(data=counts,tree=tree)
hillsim(data=counts,dist=dist)

hillsim(data=counts,q=0)
hillsim(data=counts,metric=c("C","U"),q=c(0,0.5,1),tree=tree)
hillsim(data=counts,q=2,dist=dist)
```

### Hill numbers dissimilarity metrics

The function ***hilldiss()*** computes overall dissimilarity metrics based on the Hill numbers beta diversity following Chiu et al. (2014). By default the function outputs the following four dissimilarity metrics, although users can limit the output to a single or some of those metrics. Calculations can be done for a single q-value or multiple q-values. Note these metrics are the complements (1-X) of the similarity metrics computed by ***hillsim()***. Also note that this function does not yield useful data for beta-diversity based ordinations, ***hillpair()*** should be used instead.

- **S** ***(Jaccard-type turnover):*** it quantifies the normalised species turnover rate in a sample relative to the total pool of samples (i.e., gamma diversity).
- **V** ***(Sørensen-type turnover):*** it quantifies the normalised species turnover rate in a sample relative to one sample (i.e., alpha diversity).
- **U** ***(Jaccard-type overlap-complement):*** it quantifies the proportion of non-shared species in the total pool of samples. Therefore, this metric quantifies dissimilarity from the perspective of the pool of samples.
- **C** ***(Sørensen-type overlap-complement):*** it quantifies the effective average proportion of non-shared OTUs/ASVs/MAGs in samples. Therefore, this metric quantifies dissimilarity from the perspective of a single sample.

```r
hilldiss(data=counts)
hilldiss(data=counts,tree=tree)
hilldiss(data=counts,dist=dist)

hilldiss(data=counts,q=0)
hilldiss(data=counts,metric=c("C","U"),q=c(0,0.5,1),tree=tree)
hilldiss(data=counts,q=2,dist=dist)
```

### Hill numbers pairwise dissimilarities

Using ***hillpair()***, dissimilarities are computed for every pairwise combination of samples, rather than the entire sample set (see ***hilldiss()***).

```r
hillpair(data=counts)
hillpair(data=counts,tree=tree)
hillpair(data=counts,dist=dist)

hillpair(data=counts,q=0)
hillpair(data=counts,q=1,tree=tree)
hillpair(data=counts,q=2,dist=dist)
```

#### Using pairwise dissimilarities in an ordination

Pairwise dissimilarity values of any type of diversity (neutral, phylogenetic or functional) and q-value (q=0, q=1, q=2, etc.) can be used to generate dissimilarity-based ordinations. Each type of diversity and q-value provides a different type of information, which is why it is recommendable to plot multiple ordinations to identify the main features contributing to the separation across contrasted groups.

```r
library(vegan)
library(ggplot2)

#Remove outliers to avoid distorting ordination
hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=1)

#Other Hill numbers can also be used, e.g.:
#hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=0)
#hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=1,tree=tree)
#hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=1,dist=dist)

# Generate NMDS ordination
hill_pair_dis_nmds <- hill_pair_dis %>%
				select(first,second,C) %>% #based on dissimilarity metric C
				as.data.frame() %>%
        list2dist() %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample")

#Add metadata
hill_pair_dis_nmds <- hill_pair_dis_nmds %>%
      left_join(metadata, by = join_by(sample == sample))

#Plot ordination
ggplot(hill_pair_dis_nmds, aes(x=NMDS1,y=NMDS2, color=population)) +
        geom_point(size=3) +
        scale_color_manual(values = c("#E3D97B","#46edc8","#374d7c")) +
        theme_classic() +
        theme(legend.position="bottom", legend.box="vertical")
```
The above code will result in the below NMDS ordination. Samples from population 1 differ considerably from samples in the other two populations, although some degree of separation is still observed between populations 2 and 3.
 ![NMDS ordination of neutral diversities of q=1](/images/nmds.png)

The results look slightly different when phylogenetic Hill numbers are used. Populations 2 and 3 are not separated as well as in the previous analysis, probably because the differences between both populations were driven by phylogenetically related taxa.
```r
hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=1,tree=tree)
```
 ![NMDS ordination of phylogenetic diversities of q=1](/images/nmds2.png)

 The results look very different when neutral Hill numbers of q=0 (richness) are used instead. When relative abundances are not taken into account, the population structure is lost. This example highlights the importance of using multiple combinations of Hill numbers to understand the diversity structure of samples and populations.
 ```r
hill_pair_dis <- hillpair(data=counts[,-c(12,16,19)],q=0)
 ```
  ![NMDS ordination of phylogenetic diversities of q=1](/images/nmds3.png)

### Hill numbers redundancy metrics

Using ***hillred()*** phylogenetic and functional redundancy estimations can be carried out relying on the following increasing form of an exponential decay function: y=a*2^(-x/b)+c, where:
**a** is the range of the functional Hill numbers (y value) between the lowest value in the set samples or population and the saturation value
**b** is the units of neutral Hill numbers required to increase half of the range of functional Hill numbers, and the main parameter of interest.
**c** is the value of the saturation point of functional Hill numbers (y value).

As b is dependent of the maximum neutral diversity, the redundancy metric is derived from the division between the estimated parameter b and the maximum neutral diversity value in the sample set.

The function outputs a table with the redundancy value and the estimated parameters, which can be used to plot redundancy charts.

```r
hillred(data=counts,tree=tree)
hillred(data=counts,dist=dist)

hillred(data=counts,q=1,tree=tree)
hillred(data=counts,q=2,dist=dist)
```

#### Phylogenetic redundancy

Example of phylogenetic redundancy:

```r
redundancy <- hillred(data=counts,q=1,tree=tree)
xydata <- as.data.frame(t(rbind(hilldiv(data=counts,q=1),hilldiv(data=counts,q=1,tree=tree))))
colnames(xydata) <- c("neutral","phylogenetic")

#Declare relationship function
relationship_function <- function(x, a, b, c) {return(-a*2^(-x/b)+c)}

#Plot data points
ggplot(xydata, aes(x = xydata[,1], y = xydata[,2])) +
  geom_point() +
  geom_smooth(method = "nls", formula = y ~ relationship_function(x, a, b, c), method.args = list(start = redundancy[,-1]), se = FALSE, colour="#E3D97B") +
  labs(x = "Neutral diversity", y = "Phylogenetic diversity") +
  theme_classic() +
  annotate('text', label=paste0(" Redundancy: ",round(redundancy[,1],2)), x=-Inf, y=Inf, hjust=0, vjust=1)
```

![Phylogenetic redundancy of q=1](/images/phylo_redundancy.png)

#### Functional redundancy

Example of functional redundancy:

```r
redundancy <- hillred(data=counts,q=1,dist=dist)
xydata <- as.data.frame(t(rbind(hilldiv(data=counts,q=1),hilldiv(data=counts,q=1,dist=dist))))
colnames(xydata) <- c("neutral","functional")

#Declare relationship function
relationship_function <- function(x, a, b, c) {return(-a*2^(-x/b)+c)}

#Plot data points
ggplot(xydata, aes(x = xydata[,1], y = xydata[,2])) +
  geom_point() +
  geom_smooth(method = "nls", formula = y ~ relationship_function(x, a, b, c), method.args = list(start = redundancy[,-1]), se = FALSE, colour="#46edc8") +
  labs(x = "Neutral diversity", y = "Functional diversity") +
  theme_classic() +  
  annotate('text', label=paste0(" Redundancy: ",round(redundancy[,1],2)), x=-Inf, y=Inf, hjust=0, vjust=1)

```

![Functional redundancy of q=1](/images/func_redundancy.png)


## References

- Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.
- Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.
- Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88(10), 2427–2439.
- Chao, A., Chiu, C.-H., & Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society of London. Series B, Biological Sciences, 365(1558), 3599–3609.
- Chao, A., Chiu, C.-H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93(9), 2037–2051.
- Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.
- Chao, A., Chiu, CH., Jost, L. (2016). Phylogenetic Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers. In: Pellens, R., Grandcolas, P. (eds) Biodiversity Conservation and Phylogenetic Systematics. Topics in Biodiversity and Conservation, vol 14. Springer, Cham.
- Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.
- Alberdi A., Gilbert M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19(4), 804-817.

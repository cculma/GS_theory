# This is the temp repositiory to implement genomic selection #


## covariance ##

In genetics, covariance serves a basis for computation of kinship matrix (genetic relationship matrix (GRM)) enabling inference on population structure from sample with no known close relatives as well as inference on estimation of heritability of complex tratis

## Multi-trait genomic selection (MTGS) ##
MTGS assings weights to different traits relative to their economic importance: index selection. Classical index selection only optimizes genetic gain in the next generation and requires optimizing nonlinear breeding objectives.

Formulas (1) to (3) are just to remember the mathematical expressions in $\LaTeX$ (https://www.overleaf.com/learn/latex/Subscripts_and_superscripts).

$$x^2+y^2=z^2$$(1)

$$a_1^2 + a_2^2 = a_3^2$$(2)

$$\sum_{i=1}^{\infty} \frac{1}{n^s} = \prod_p \frac{1}{1 - p^{-s}}$$(3)

$$ y=\sum_{i}g_ix_i+\sum_{ij}g_ig_jz_{ij}+O(g^3)+noise$$(4)

(4) is the general model for quantitative phenotype where:
$y$ = inividual phenotype

$g_i$ = individual genotype (e.g. list of SNPs)

$x_i$ = linear effect sizes

$z_{ij}$ = tensors of non linear effect sizes 

Approach maximizing certain traits while keeping others within desirable ranges.
LAS = look-ahead selection algorithm
LAMS = look-ahead mate selection 
MTLAS is more effective balancing multiple traits compared with index selection
GEBV = genomic estimated breeding values are calculated as the sum of the estimated marker effects
Goal 1
Development of models to improve the accuracy of GEBV prediction
These studies have been focused on selecting the parents of the next generation by defining new quantitative selection merits
Mate selection optimize the contribution of potential parents to the next generation based on maximizing a desired breeding objective
The optimization is performed with respect to the next generation
LAMS optimize parental contributions with respect to grand-progeny (two generations in the future)
ST-GS = single trait genomic selection
Pareto efficiency is a state where resources are allocated in the most efficient manner
The accuracy of genomic prediction critically depends on the level of relatedness between the training population and the testing population
The use of multivariate linear mixed models in genetic evaluations provides a basis for inference about traits' integration
 
## PLS-CA ##

GS using all markers
GS using marker selection
Use of BLUEs estimates are used with genomic data to estimate marker weights prior to cosntruction of the marker based relationship matrix

PLS-CA: Partial Least Squares Canonical Analysis

<img src="https://render.githubusercontent.com/render/math?math=\sum_{i=1}^{\infty} \frac{1}{n^s} = \prod_p \frac{1}{1 - p^{-s}}">(3)


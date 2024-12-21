# BayesTraj
BayesTraj is a probabilistic pseudotime framework. BayesTraj

- **A.** Inputs to BayesTraj

- **B.** Illustration of the Bayesian Model

- **C.** Trend of Marker Gene Changes During Differentiation
  
- **D.** Generation of Pseudo-time and Branch Probabilities

- **E.** Branch-specific gene detection

  ![fig1_画板 1](https://github.com/user-attachments/assets/557a698d-6c90-410a-aaf8-1045d956a93c)



# Getting started

## Installation

```r
# install.packages("devtools")
devtools::install_github("SDU-W-ZhangLab/BayesTraj")
```

## Model fitting

### linear trajectories

To evaluate BayesTraj's ability to infer cellular trajectories, we applied it to time-series scRNA-seq datasets of **mESC and hESC differentiation into endoderm cells**. These datasets included multiple time points with temporal annotations, allowing us to assess the inferred pseudotime. We tested whether BayesTraj, **using small panels of marker genes**, could accurately reconstruct the differentiation trajectory. 

<p align="center">—— <strong>mesc</strong> ——</p>

(1) 筛选基因
(2) 输入marker基因表达和分支结构+参数介绍
(3) 输出结果+输出框架结构介绍

<p align="center">—— <strong>hesc</strong> ——</p>


### bifurcation trajectories

We use scRNA-seq datasets with bifurcation structures to examine whether BayesTraj can successfully detect branching trajectories. One such dataset is **a mouse lung dataset** containing 185 cells collected at E14.5, E16.5, and E18.5, which differentiate into type I (AT1) and type II (AT2) pneumocytes.



### tree trajectories

We utilized a scRNA-seq dataset from **the mouse hematopoietic system**, encompassing all developmental stages, to further evaluate BayesTraj's ability to reconstruct multiple lineages.


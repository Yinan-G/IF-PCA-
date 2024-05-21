# IF-PCA Plus

This repository contains the code used in the manuscript Comparison of Clustering Methods on Single-Cell
Data.

## Datasets
Two example datasets are provided, one for Microarray and one for single-cell data.

The 8 single cell datasets can be download at https://data.mendeley.com/drafts/nv2x6kf5rd and the 10 microarray data sets can be downloaded at https://data.mendeley.com/datasets/cdsz2ddv3t.

## Codes
The implementation of IF-PCA+ consists of main two parts: a manifold fitting part and an IF-PCA part.

The manifold fitting method is implemented in Matlab. The original code is from Professor Zhigang Yao's scAMF paper.
The IF-PCA part is implemented in Python.

To run IF-PCA+ on a selected dataset:
1. Open the "Yao2_denoising.m" file from Matlab and select the desired dataset by modifying the name varaible.
If the dataset selected contains negative entries, comment out the codes corresponding to log transformation. The denoised data array will be saved under "datasetname_Y2_denoised_data.mat"

2. Open the IFPCAplus file in Python and run the code accordingly.





External links:
The code to implement IF-PCA, IF-VAE, Seurat, SC3 can be find at https://github.com/ZhengTracyKe/IFPCA. The code to implement scAMF can be find at https://github.com/zhigang-yao/scAMF. 



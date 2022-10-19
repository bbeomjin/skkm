# skkm
Sparse Kernel k-Means Clustering

```skkm``` is an R package. ```skkm``` provides functions for fitting the sparse kernel k-means clustering.

## 1. INSTALLATION

(1) From GitHub
```{r}
> library(devtools)
> install_github("bbeomjin/skkm")
```

## 2. USAGE NOTES

(1) Description of R functions in ```skkm```

- Descriptions of arguments in the functions in ```skkm``` can be obtained by help() or ? in R prompt, and documentation of ```skkm```.   

(2) List of R functions in ```skkm``` package

- ```tune.skkm``` : A function for tuning the hyper-parameter of the sparse kernel k-means clustering.

- ```skkm``` : A function to execute the sparse kernel k-means clustering.

- ```tune.kkmeans``` : A function for tuuing the kernel parameter of kernel k-means clustering.

- ```kkmeans``` : A function to perform the kernel k-means clustering.

## 3. Examples
```{r}
require(fossil)

# Generate the simulated data through the Smiley scenario
n = 100
p = 2
dat = generateSmiley(n = n, p = p, seed = 1, with_noise = TRUE, noise_p = 5)

# Grid for the Gaussian kernel parameter
sigma = c(0.25, 0.5, 0.75, 1)

# Tuning parameters for the Sparse kernel k-means clustering with the exact search
# The nCores argument only works on Linux
# For Windows, set nCores = 1 
tuned_skkm = tune.skkm(x = dat$x, nCluster = nclusters, s = NULL, ns = 20, nPerms = 25,
                       nStart = 1, kernel = "gaussian-2way", kparam = sigma, search = "exact", 
                       opt = TRUE, nInit = 20, nCores = 20)

# Tuning parameters for the Sparse kernel k-means clustering with the binary search
# Does not excecute
# tuned_skkm = tune.skkm(x = dat$x, nCluster = nclusters, s = NULL, ns = 20, nPerms = 25,
#                        nStart = 1, kernel = "gaussian-2way", kparam = sigma, search = "binary", 
#                        opt = TRUE, nInit = 20, nCores = 20)

skkm_clusters = tuned_skkm$optModel$optClusters
ari_skkm = adj.rand.index(dat$y, skkm_clusters)
```

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

- ```skkm``` : A function to execute the sparse kernel k-means algorithm.

# ECLRMC

Ensemble correlation-based low-rank matrix completion method (ECLRMC) is an extension to the LRMC based methods. Traditionally, the LRMC based methods give identical importance to the whole data which results in emphasizing on the commonality of the data and overlooking the subtle but crucial differences. 
This method aims to overcome the equality assumption problem that exists in the current LRMS based methods. Ensemble correlation-based low-rank matrix completion (ECLRMC) takes consideration of the specific characteristic of each sample and performs LRMC on the set of samples with a strong correlation. 
It uses an ensemble learning method to improve the imputation performance. Since each sample is analyzed independently this method can be parallelized by distributing imputation across many computation units or GPU platforms. This package provides three different methods (LRMC, CLRMC and ECLRMC) for data imputation. 
There is also an NRMS function for evaluating the result.
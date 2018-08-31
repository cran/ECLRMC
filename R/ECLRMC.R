#' Takes an incomplete matrix and returns the imputed matrix using LRMC method.
#' @param x An m by n matrix with NAs
#' @return An m by n matrix with imputed values
#' @references Chen, Xiaobo, et al. "Ensemble correlation-based low-rank matrix completion with applications to traffic data imputation." Knowledge-Based Systems 132 (2017): 249-262.
#' @examples
#' x = matrix(c(5.1, 4.9, NA, 4.6, 3.5, 3.0, 3.2, 3.1, 1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
#' LRMC(x)
#' @export
LRMC <- function(x) {
    #----------------------------------------------------------------
    # Sample x
    #----------------------------------------------------------------
    InputMatrixRows = nrow(x)
    InputMatrixColumns = ncol(x)
    
    #----------------------------------------------------------------
    # LRMC
    #----------------------------------------------------------------
    mySvd = softImpute(x = x , type = "svd")
    lrmcMatrix = complete(x, mySvd)
    return(lrmcMatrix)
}

#' Takes an incomplete matrix and returns the imputed matrix using CLRMC method.
#' @param x An m by n matrix with NAs
#' @param beta A value in [0,1] range. Higher beta value means comparing each row with more nearest neighbours. Default value = 0.1
#' @return An m by n matrix with imputed values
#' @references Chen, Xiaobo, et al. "Ensemble correlation-based low-rank matrix completion with applications to traffic data imputation." Knowledge-Based Systems 132 (2017): 249-262.
#' @examples
#' x = matrix(c(5.1, 4.9, NA, 4.6, 3.5, 3.0, 3.2, 3.1, 1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
#' CLRMC(x, beta = 0.2)
#' @export
CLRMC <- function(x, beta = 0.1) {
  #----------------------------------------------------------------
  # Sample x
  #----------------------------------------------------------------
  # x = matrix(c(5.1, 4.9, NA, 4.6, 
  #              3.5, 3.0, 3.2, 3.1, 
  #              1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
  InputMatrixRows = nrow(x)
  InputMatrixColumns = ncol(x)
  
  #----------------------------------------------------------------
  # Step 1 - LRMC
  #----------------------------------------------------------------
  
  lrmcMatrix = LRMC(x)
  
  #----------------------------------------------------------------
  # Step 2 - Distance Matrix
  #----------------------------------------------------------------
  distanceMatrix = data.frame()
  
  # Weighted Matrix
  weights = array(
    data = NA,
    dim = c(InputMatrixRows, InputMatrixRows, InputMatrixColumns)
  )
  for (i in  1:InputMatrixRows)
  {
    for (j in  1:InputMatrixRows)
    {
      for (n in  1:InputMatrixColumns)
      {
        if (is.na(x[i, n]) || is.na(x[j, n]))
          weights[i, j, n] = 0.1
        else
          weights[i, j, n] = 1
      }
    }
  }
  
  for (i in  1:InputMatrixRows)
  {
    for (j in  1:InputMatrixRows)
    {
      weightSum = sum(weights[i, j, ])
      weightNormalized = c()
      XiNormalized = c()
      XjNormalized = c()
      SumWeightedXi = 0
      SumWeightedXj = 0
      for (n in  1:InputMatrixColumns)
      {
        weightNormalized[n] = weights[i, j, n] / weightSum
        SumWeightedXi = SumWeightedXi + (weightNormalized[n] * lrmcMatrix[i, n])
        SumWeightedXj = SumWeightedXj + (weightNormalized[n] * lrmcMatrix[j, n])
      }
      for (n in  1:InputMatrixColumns)
      {
        XiNormalized[n] = lrmcMatrix[i, n] - SumWeightedXi
        XjNormalized[n] = lrmcMatrix[j, n] - SumWeightedXj
      }
      
      v1 = sum(weightNormalized[] * XiNormalized[] * XjNormalized[])
      v2 = sum(weightNormalized[] * (XiNormalized[] ^ 2))
      v3 = sum(weightNormalized[] * (XjNormalized[] ^ 2))
      
      if (v1==0)
        cor = 0
      else
        cor = v1 / sqrt(v2 * v3)
      
      absValue = round(abs(cor), 8)
      distanceMatrix[i, j] = 1 - absValue
    }
  }
  
  #----------------------------------------------------------------
  # Step 3 - Adaptive KNN
  #----------------------------------------------------------------
  # Input Parameter

  orderedDistanceMatrix = t(apply(distanceMatrix, 1, sort))
  distanceRatioMatrix = data.frame(matrix(NA, InputMatrixRows, InputMatrixRows))
  distanceNeighborsNumber = c()
  distanceNeighborsValue = c()
  for (i in  1:InputMatrixRows)
  {
    distanceNeighborsNumber[i] = 0
    distanceNeighborsValue[i] = 0
    booldistanceNeighbors = FALSE
    
    v2 = sum(orderedDistanceMatrix[i, ])
    
    for (p in  1:InputMatrixRows)
    {
      v1 = sum(orderedDistanceMatrix[i, 1:p])
      distanceRatioMatrix[i, p] = v1 / v2
      
      if (distanceRatioMatrix[i, p] >= beta &&
          booldistanceNeighbors == FALSE)
      {
        distanceNeighborsNumber[i] = p
        distanceNeighborsValue[i] = distanceRatioMatrix[i, p]
        booldistanceNeighbors = TRUE
      }
    }
  }
  
  #----------------------------------------------------------------
  # Step 4 - Correlation-Based LRMC
  #----------------------------------------------------------------
  knnLRMC = list(matrix(NA, InputMatrixRows, InputMatrixColumns),
                 InputMatrixRows)
  knnIndex = data.frame(matrix(NA, nrow = InputMatrixRows + 1 , ncol = InputMatrixColumns))
  
  for (i in 1:InputMatrixRows) {
    kNN = data.frame(matrix(NA, nrow = distanceNeighborsNumber + 1 , ncol = InputMatrixColumns))
    
    for (j in 1:InputMatrixColumns) {
      kNN[1, j] = x[i, j]
    }
    
    distanceCoordinate <-
      c(apply(distanceMatrix[i, ], 1, function(x)
        order(x)[1:distanceNeighborsNumber[i]]))
    
    for (k in 1:distanceNeighborsNumber[i]) {
      for (j in 1:InputMatrixColumns) {
        kNN[k, j] = x[distanceCoordinate[k], j]
      }
      knnIndex[k, i] = distanceCoordinate[k]
    }
    
    mySvd2 = softImpute(x = kNN , type = "svd")
    knnLRMC[[i]] = complete(kNN, mySvd2)
  }
  
  clrmc = c()
  for (i in 1:InputMatrixRows)
  {
    clrmc = rbind(clrmc, knnLRMC[[i]][1,])
  }
  clrmc = data.matrix(clrmc)
  return(clrmc)
}

#' Takes an incomplete matrix and returns the imputed matrix using ECLRMC method.
#' @param x An m by n matrix with NAs
#' @param beta A value in [0,1] range. Higher beta value means comparing each row with more nearest neighbours. Default value = 0.1
#' @return An m by n matrix with imputed values
#' @references Chen, Xiaobo, et al. "Ensemble correlation-based low-rank matrix completion with applications to traffic data imputation." Knowledge-Based Systems 132 (2017): 249-262.
#' @examples
#' x = matrix(c(5.1, 4.9, NA, 4.6, 3.5, 3.0, 3.2, 3.1, 1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
#' ECLRMC(x, beta = 0.2)
#' @export
ECLRMC <- function(x, beta = 0.1){
  #----------------------------------------------------------------
  # Sample x
  #----------------------------------------------------------------
  # x = matrix(c(5.1, 4.9, NA, 4.6, 
  #              3.5, 3.0, 3.2, 3.1, 
  #              1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
  InputMatrixRows = nrow(x)
  InputMatrixColumns = ncol(x)
  
  #----------------------------------------------------------------
  # Step 1 - LRMC
  #----------------------------------------------------------------
  
  lrmcMatrix = LRMC(x)
  
  #----------------------------------------------------------------
  # Step 2 - Distance Matrix
  #----------------------------------------------------------------
  distanceMatrix = data.frame()
  
  weights = array(
    data = NA,
    dim = c(InputMatrixRows, InputMatrixRows, InputMatrixColumns)
  )
  for (i in  1:InputMatrixRows)
  {
    for (j in  1:InputMatrixRows)
    {
      for (n in  1:InputMatrixColumns)
      {
        if (is.na(x[i, n]) || is.na(x[j, n]))
          weights[i, j, n] = 0.1
        else
          weights[i, j, n] = 1
      }
    }
  }
  
  for (i in  1:InputMatrixRows)
  {
    for (j in  1:InputMatrixRows)
    {
      weightSum = sum(weights[i, j, ])
      weightNormalized = c()
      XiNormalized = c()
      XjNormalized = c()
      SumWeightedXi = 0
      SumWeightedXj = 0
      for (n in  1:InputMatrixColumns)
      {
        weightNormalized[n] = weights[i, j, n] / weightSum
        SumWeightedXi = SumWeightedXi + (weightNormalized[n] * lrmcMatrix[i, n])
        SumWeightedXj = SumWeightedXj + (weightNormalized[n] * lrmcMatrix[j, n])
      }
      for (n in  1:InputMatrixColumns)
      {
        XiNormalized[n] = lrmcMatrix[i, n] - SumWeightedXi
        XjNormalized[n] = lrmcMatrix[j, n] - SumWeightedXj
      }
      v1 = sum(weightNormalized[] * XiNormalized[] * XjNormalized[])
      v2 = sum(weightNormalized[] * (XiNormalized[] ^ 2))
      v3 = sum(weightNormalized[] * (XjNormalized[] ^ 2))
      if (v1==0)
        cor = 0
      else
        cor = v1 / sqrt(v2 * v3)
      absValue = round(abs(cor), 8)
      distanceMatrix[i, j] = 1 - absValue
    }
  }
  
  #----------------------------------------------------------------
  #Step 3 - Adaptive KNN
  #----------------------------------------------------------------
  orderedDistanceMatrix = t(apply(distanceMatrix, 1, sort))
  distanceRatioMatrix = data.frame(matrix(NA, InputMatrixRows, InputMatrixRows))
  distanceNeighborsNumber = c()
  distanceNeighborsValue = c()
  for (i in  1:InputMatrixRows)
  {
    distanceNeighborsNumber[i] = 0
    distanceNeighborsValue[i] = 0
    booldistanceNeighbors = FALSE
    
    v2 = sum(orderedDistanceMatrix[i, ])
    
    for (p in  1:InputMatrixRows)
    {
      v1 = sum(orderedDistanceMatrix[i, 1:p])
      distanceRatioMatrix[i, p] = v1 / v2
      if (distanceRatioMatrix[i, p] >= beta &&
          booldistanceNeighbors == FALSE)
      {
        distanceNeighborsNumber[i] = p
        distanceNeighborsValue[i] = distanceRatioMatrix[i, p]
        booldistanceNeighbors = TRUE
      }
    }
  }
  
  #----------------------------------------------------------------
  #Step 4 - Correlation-Based LRMC
  #----------------------------------------------------------------
  knnLRMC = list(matrix(NA, InputMatrixRows, InputMatrixColumns),
                 InputMatrixRows)
  # knnIndex
  knnIndex = data.frame(matrix(NA, nrow = InputMatrixRows + 1 , ncol = InputMatrixColumns))
  for (i in 1:InputMatrixRows) {
    kNN = data.frame(matrix(NA, nrow = distanceNeighborsNumber + 1 , ncol = InputMatrixColumns))
    for (j in 1:InputMatrixColumns) {
      kNN[1, j] = x[i, j]
    }
    distanceCoordinate <-c(apply(distanceMatrix[i, ], 1, function(x) order(x)[1:distanceNeighborsNumber[i]]))
    for (k in 1:distanceNeighborsNumber[i]) {
      for (j in 1:InputMatrixColumns) {
        kNN[k, j] = x[distanceCoordinate[k], j]
      }
      knnIndex[k, i] = distanceCoordinate[k]
    }
    
    mySvd2 = softImpute(x = kNN , type = "svd")
    knnLRMC[[i]] = complete(kNN, mySvd2)
  }
  
  #----------------------------------------------------------------
  #Step 5 - Ensemble Learning
  #----------------------------------------------------------------
  
  standardDeviation = c()
  finalMatrix = matrix(0, InputMatrixRows, InputMatrixColumns)
  duplicateRowsMatrix = matrix(0, InputMatrixRows, InputMatrixRows)
  duplicateRowsNumber = c()
  duplicateRowsDistanceMatrix = matrix(0, InputMatrixRows, InputMatrixRows)
  
  for (i in 1:InputMatrixRows) {
    counter = 1
    duplicateRowsDistances = c()
    for (j in 1:InputMatrixRows) {
      existRow = which(knnIndex[, j] == i)[1]
      if (!is.na(existRow))
      {
        duplicateRowsMatrix[i, counter] = j
        duplicateRowsDistances = c(duplicateRowsDistances, distanceMatrix[i, j])
        counter = counter + 1
      }
    }
    
    a = sd(duplicateRowsDistances)
    if(is.na(a)) {
      standardDeviation[i] = 0
    } else {
      standardDeviation[i] = sd(duplicateRowsDistances)
    }
    s = standardDeviation[i]
    duplicateRowsNumber[i] = counter - 1
  }
  
  rWeight = matrix(0, InputMatrixRows, InputMatrixRows)
  for (i in 1:InputMatrixRows) {
    for (p in 1:duplicateRowsNumber[i])
    {
      distanceIP = distanceMatrix[i, duplicateRowsMatrix[i, p]]
      if (standardDeviation[i]==0)
      {
        rWeight[i,p] = 1
      }
      else
      {
        rWeight[i,p] = exp(-((distanceIP ^ 2) / (standardDeviation[i]^2)))
      }
    }
  }
  
  for (i in 1:InputMatrixRows) {
    for (j in 1:InputMatrixColumns)
    {
      s = standardDeviation[i]
      for (p in 1:duplicateRowsNumber[i])
      {
        knnNumber = strtoi(duplicateRowsMatrix[i, p])
        KnnRow = strtoi(which(knnIndex[, knnNumber] == i))
        rWeightValue = rWeight[i,p] / sum(rWeight[i,])
        if (length(knnNumber) == 0){
          finalMatrix[i, j] = lrmcMatrix[i,j]
        } else if (knnNumber==0)
        {
          finalMatrix[i, j] = lrmcMatrix[i,j]
        }
        else
        {
          finalMatrix[i, j] = finalMatrix[i, j] + rWeightValue * knnLRMC[[knnNumber]][KnnRow, j]
        }
      }
    }
  }
  return(finalMatrix)
}

#' Normalized Root Mean Square (NRMS) value of two matrices for evaluating their similarity (lower is better)
#' @param imputed An m by n matrix
#' @param original An m by n matrix
#' @return Returns the NRMS value of the given matrices
#' @examples
#' x = matrix(c(5.1, 4.9, NA, 4.6, 3.5, 3.0, 3.2, 3.1, 1.4, NA, 1.3, 1.5), byrow = TRUE, ncol=4)
#' a = ECLRMC(x, beta = 0.2)
#' b = LRMC(x)
#' NRMS(a,b)
#' @export
NRMS <- function(imputed, original){
  numerator = sqrt(sum((imputed-original)^2))
  if(numerator == 0) {
    return(0)
  }
  result = numerator/sqrt(sum(original^2))
  return(result)
}
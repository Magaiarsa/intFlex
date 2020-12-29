## condense an occupancy matrix into as few columns as possible
condense.mat <- function(m, yrs) {
  condense.year <- function(y) {
    condense <- function(x, y) {
      tmp <- rep(NA, length(yrs==y))
      x <- x[yrs==y]
      if(any(!is.na(x)))
        tmp[1:sum(!is.na(x))] <- x[!is.na(x)]
      tmp
    }
    m.y <- t(apply(m, 1, condense, y=y))
    keep <- apply(m.y, 2, function(x) any(!is.na(x)))
    if(sum(keep)>0) {
      m.y <- matrix(m.y[,keep], ncol=sum(keep))
      colnames(m.y) <- rep(y, ncol(m.y))
      rownames(m.y) <- rownames(m)
    }
    if(sum(keep)==0) {
      m.y <- m.y[,keep]
    }
    m.y
  }
  by.year <- lapply(unique(yrs), condense.year)
  do.call(cbind, by.year)
}

## make 4D matrix (site x year x replicate x species)
make.mat <- function(mat, threshold, nzero) {
  
  ## set all >1 entries to 1
  mat[mat>1] <- 1

  ## drop species observed less than 'threshold' times in total
  mat <- mat[,,apply(mat, 3, sum, na.rm=TRUE)>=threshold]
  
  ## drop columns that are NA for all species
  mat <- mat[,apply(mat, 2, function(x) any(!is.na(x))),]
  nms <- dimnames(mat)

  ## augment data (e.g., add matrices of missing species)
  if(nzero>0) {
    n.entries <- apply(mat, 3, function(x) sum(x>=0, na.rm=TRUE))
    ii <- which(n.entries==max(n.entries))[1]
    mat.zero <- array(mat[,,ii]*0, dim=c(c(dim(mat)[1:2], nzero)))
    nms$species <- c(nms$species, sprintf("sp.%d", 1:nzero))
    mat <- array(c(mat, mat.zero),
                 dim=c(dim(mat)[1:2],dim(mat)[3]+nzero),
                 dimnames=nms)
  }

  Y <- as.numeric(sapply(strsplit(dimnames(mat)$date, split="-"),
                         function(x) x[[1]]))
  
  get.date.mat <- function(sp) {
    dates <-
      lapply(1:nrow(mat[,,sp]), function(x)
             colnames(mat[,,sp])[!is.na(mat[x,,sp])])
    names(dates) <- dimnames(mat)$site
    dates
  }
  dates <- lapply(1:(dim(mat)[3]), get.date.mat)
  names(dates) <- dimnames(mat)$species

  nms$date <- colnames(condense.mat(mat[,,1], yrs=Y))
  arr.tmp <- apply(mat, 3, condense.mat, yrs=Y)

  if(is.list(arr.tmp)) {
    max.dim <- max(sapply(arr.tmp, dim)[2,])
    augment <- function(m) {
      if(dim(m)[2]==max.dim) return(m)
      cbind(m, matrix(NA, nrow=nrow(m), ncol=max.dim-ncol(m)))
    }
    arr.tmp <- sapply(arr.tmp, augment)
  }

  mat <- array(arr.tmp, dim=c(dim(mat)[1],
                          dim(arr.tmp)[1]/dim(mat)[1],
                          dim(mat)[3]),
               dimnames=nms)
  list(mat=mat, dates=dates)
}

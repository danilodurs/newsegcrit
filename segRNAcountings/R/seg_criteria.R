#' log-transformation function for segmentation
#'
#' @param rna.data : a vector of counts from RNA-sequencing
#'
#' @return a vector of log-transformed data
#'
#'
#' @examples
#' log.data <- log.transform(dataset1)
#' plot(dataset1, type="l")
#' plot(log.data, type="l")
log.transform <- function(rna.data)
{
  log.data <- log(rna.data+1)
  n <- length(rna.data)
  ##transformation "HALL"
  x   <- log.data
  wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
  mat <- wei %*% t(x)
  mat[2, -n] = mat[2, -1]
  mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
  mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]
  est.sd <- sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
  return(log.data/est.sd)
}


#' returning the index of a penalty among a vector of intermediate penalties
#'
#' @param lambda a penalty value
#' @param o.crops the output object from CROPS.RFPOP_
#'
#' @return the index of the penalty, an integer
#'
#' @examples
#' log.t1 <- log.transform(dataRNA[,1])
#' crops.out  <- CROPS.RFPOP_(log.t1,min_pen = 5,max_pen = 10,lthreshold = 3)
#' getInd <- getIndex(lambda = 10*log(n), o.crops = crops.out)
#' getInd
#'
getIndexLambda <- function(lambda, o.crops){
  n.seg <- ncol(o.crops[[1]])
  index <- sum( o.crops[[1]][2, ] <= lambda)
  if(index == 0)     return(1)
  return(index)
}

#' Segmenting a list of dataset for a given segmentation algo and a penalty range
#'
#' @param list.rna.data a list of dataset
#' @param penalty_range a vector of length 2. Respectively. The min and max penalties
#' @param mc.cores a paramater for mclapply, 1 by default... (CF the appropriate documentation about mcapply...)
#' @return a list of values: a list of tau for each beta, a list of smt for each dataset, a list of intermediate penalties.
#'
#' @examples
#' l.d1     <- log.transform(dataset1)
#' l.d2     <- log.transform(dataset2)
#' l.d3     <- log.transform(dataset3)
#' l.data   <- list(l.d1,l.d2,l.d3)
#' list.res <- seg.func(list.rna.data=l.data,penalty_range=c(10,100),mc.cores = 1)
#' View  (list.res)
#'
seg.func <- function(list.rna.data, penalty_range, mc.cores=1)
{
  #First, crops
  crops <- mclapply(list.rna.data, FUN=function(rna.data) {
    CROPS.RFPOP_(data = rna.data , min_pen = penalty_range[1] , max_pen = penalty_range[2])
  }, mc.cores=mc.cores)
  ## getAllPenalties
  all.penalties <- unlist(lapply(crops, FUN=function(x) x[[1]][2, ]))
  intermediate.penalties <- sort(unique(all.penalties))

  list.tau.results <- lapply(crops, FUN=function(o.crops){
    res <- lapply(intermediate.penalties, FUN=function(lambda){
      i <- getIndexLambda(lambda, o.crops)
      o.crops[[2]][[i]]
    })
  })

  list.smt.results <- lapply(crops, FUN=function(o.crops){
    res <- lapply(intermediate.penalties, FUN=function(lambda){
      i <- getIndexLambda(lambda, o.crops)
      o.crops[[3]][[i]]
    })
  })

  list.results <- list(list.tau.results , list.smt.results, intermediate.penalties)
  class(list.results) <- "SEGMENTATIONS"
  return(list.results)
}

#' Computing MSE depending on penalty
#'
#' @param list.seg a "SEGMENTATIONS" object
#'
#' @return a matrix of results: the lines represent the observations (the replicas), the "avant-dernière" line: the mean of them. The last line: the intermediate penalties. The columns represent the mse value for each penalty.
#'
#' @examples
#' l.d1     <- log.transform(dataset1)
#' l.d2     <- log.transform(dataset2)
#' l.d3     <- log.transform(dataset3)
#' l.data   <- list(l.d1,l.d2,l.d3)
#' list.res <- seg.func(func="rob_seg",list.rna.data=l.data,penalty_range=c(10,100))
#' mse.res  <- mse.penalties(list.res)
#'
mse.penalties <- function(list.seg)
{
  list.smt.results <- list.seg[[2]]
  datasets.index   <- 1:length(list.smt.results)
  combins          <- combn(datasets.index,2,simplify = F)
  matrix.results   <- NULL
  #picking a datset of smt from each replica, for each beta
  i <- 1
  #for each beta
  while( i <= length(list.smt.results[[1]]) )
  {
    j <- 1
    datasets.list <- list()
    #defining some datasets for a given beta
    while(j <= length(list.smt.results) )
    {
      datasets.list[[j]] <- list.smt.results[[j]][[i]]
      j                  <- j+1
      cat(".")
    }
    col.mse <- NULL
    for(combin in combins)
    {
      mse <- mean((datasets.list[[combin[1]]]-datasets.list[[combin[2]]])^2)
      col.mse <- c(col.mse, mse)
      cat("+")
    }
    matrix.results <- cbind(matrix.results,col.mse)
    i <- i+1
    cat("*")
  }
  matrix.results <- rbind( matrix.results , colMeans(matrix.results))
  matrix.results <- rbind( matrix.results , list.seg[[3]])
  class(matrix.results) <- "MSE"
  return(matrix.results)
}

#' Computing NID depending on penalty
#'
#' @param list.seg an "SEGMENTATIONS" object
#'
#' @return a matrix of results: the lines represent the observations (the replicas), the "avant-dernière" line: the mean of them. The last line: the intermediate penalties. The columns represent the nid value for each penalty.
#'
#' @examples
#' l.d1     <- log.transform(dataset1)
#' l.d2     <- log.transform(dataset2)
#' l.d3     <- log.transform(dataset3)
#' l.data   <- list(l.d1,l.d2,l.d3)
#' list.res <- seg.func(func="rob_seg",list.rna.data=l.data,penalty_range=c(10,100))
#' nid.res  <- nid.penalties(list.res)
#'
nid.penalties <- function(list.seg)
{
  list.tau.results <- list.seg[[1]]
  datasets.index   <- 1:length(list.tau.results)
  combins          <- combn(datasets.index,2,simplify = F)
  matrix.results   <- NULL
  #picking a datset of smt from each replica, for each beta
  i <- 1
  #for each beta
  while( i <= length(list.tau.results[[1]]) )
  {
    j <- 1
    datasets.list <- list()
    #defining some datasets for a given beta
    while(j <= length(list.tau.results) )
    {
      datasets.list[[j]] <- list.tau.results[[j]][[i]]
      j                  <- j+1
      cat(".")
    }
    #preparing the class for each segment
    l <- 1
    all.seg.class <- list()
    for(tauset in datasets.list)
    {
      tauset <- c(0,tauset)
      #defining a segment and attribute the class
      k <- 1
      seg.class <- NULL
      while(k <= (length(tauset)-1))
      {
        seg.class <- c(seg.class,rep(k,length((tauset[k]+1):tauset[k+1])))
        k         <- k + 1
      }
      all.seg.class[[l]] <- seg.class
      l <- l+1
    }
    #####################################
    col.nid <- NULL
    for(combin in combins)
    {
      nid.res <- NID(all.seg.class[[combin[1]]],all.seg.class[[combin[2]]])
      col.nid <- c(col.nid, nid.res)
      cat("+")
    }
    matrix.results <- cbind(matrix.results,col.nid)
    i <- i+1
    cat("*")
  }
  matrix.results <- rbind( matrix.results , colMeans(matrix.results))
  matrix.results <- rbind( matrix.results , list.seg[[3]])
  class(matrix.results) <- "NID"
  return(matrix.results)
}

#' Ploting the mse depending on the penalties
#'
#' @param mse.res a "MSE" object
#'
#' @return a graph...
#' @examples plot(mse.res)
plot.MSE <- function(mse.res)
{
  plot(x=mse.res[length(mse.res[,1]),],y=mse.res[length(mse.res[,1])-1,], xlab="Penalty", ylab="MSE", type="s", col="red")
}

#' Ploting the nid depending on the penalties
#'
#' @param mse.res a "NID" object
#'
#' @return a graph...
#' @examples plot(nid.res)
#'
plot.NID <- function(nid.res)
{
  plot(x=nid.res[length(nid.res[,1]),],y=nid.res[length(nid.res[,1])-1,], xlab="Penalty", ylab="NID", type="s", col="red")
}

#' Computing different criterion for RNAs segmentations
#'
#' @param list.seg A "SEGMENTATIONS" object
#' @param criterion A string indicating the selected criterion: "MSE" or "NID" so far
#'
#' @return the selected criterion and a plot
#'
#' @examples
#' l.d1     <- log.transform(dataset1)
#' l.d2     <- log.transform(dataset2)
#' l.d3     <- log.transform(dataset3)
#' l.data   <- list(l.d1,l.d2,l.d3)
#' list.res <- seg.func(func="rob_seg",list.rna.data=l.data,penalty_range=c(10,100))
#' crit.res <- seg.criteria(list.res, criterion="MSE")
#' crit.res2 <- seg.criteria(list.res, criterion="NID")
#'
seg.criteria <- function(list.seg, criterion)
{
  crit.res <- NULL
  if(criterion == "MSE")
  {
    crit.res <- mse.penalties(list.seg)
    plot(crit.res)
  }
  else if(criterion == "NID")
  {
    crit.res <- nid.penalties(list.seg)
    plot(crit.res)
  }
  return(crit.res)
}

#' The average countings per segments for each replica
#'
#' @param list.tau list of changepoints for each replica, for a given penalty value.
#' @param list.data list of dataset for each replica.
#'
#' @return a matrix containing the average count for each segment (column), for each replica(row).
#'
#' @examples
#' l.d1 <- log.transform(dataset1)
#' seg_rob1      <- Rob_seg.std(x = l.d1,  loss = "Outlier", lambda = 25*log(length(l.d1)), lthreshold=3)
#' tau1           <- seg_rob1$t.est
#' l.d2 <- log.transform(dataset2)
#' seg_rob2      <- Rob_seg.std(x = l.d2,  loss = "Outlier", lambda = 25*log(length(l.d1)), lthreshold=3)
#' tau2           <- seg_rob2$t.est
#' l.d3 <- log.transform(dataset3)
#' seg_rob3      <- Rob_seg.std(x = l.d3,  loss = "Outlier", lambda = 25*log(length(l.d1)), lthreshold=3)
#' tau3           <- seg_rob3$t.est
#' l.data <- list(dataset1,dataset2,dataset3)
#' l.tau <- list(tau1,tau2,tau3)
#' cps <- counts.per.seg(list.tau=l.tau,list.data=l.data)
counts.per.seg <- function(list.tau,list.data)
{
  vec.tau          <- unlist(list.tau)
  vec.tau          <- sort(vec.tau[!duplicated(vec.tau)])
  vec.tau          <- c(0,vec.tau)
  mat.res          <- NULL
  for(dataset in list.data)
  {
    i            <- 1
    row.seg.mean <- NULL
    while(i < length(vec.tau))
    {
      seg          <- dataset[(vec.tau[i]+1):vec.tau[i+1]]
      seg.mean     <- mean(seg)
      row.seg.mean <- c(row.seg.mean, seg.mean)
      i            <- i+1
    }
    mat.res <- rbind(mat.res,row.seg.mean)
  }
  return(mat.res)
}

#' Segmenting a dataset using rob_seg, over a penalty range scanned by a CROPS algorithm
#'
#' @param data The dataset, a vector
#' @param min_pen minimum value of the penalty range
#' @param max_pen maximum value of the penalty range
#' @param lthreshold the threshold used to detect and rescale the outliers among the dataset (3 by default)
#'
#' @return A list: respectively a matrix of penalties (the 2nd line contains the intermediate penalties),
#'  a list of segmentations for each intermediate lambda and
#'  a list of smt for each intermediate lambda (they are obviously all the same...)
#'
#' @examples
#' log.t1 <- log.transform(dataRNA[,1])
#' crops.out  <- CROPS.RFPOP_(log.t1,min_pen = 5,max_pen = 10,lthreshold = 3)
#' View(crops.out)
#'
CROPS.RFPOP_ <- function(data, min_pen = 5, max_pen = 20, lthreshold = 3) {

  NCALC=0
  pen_interval <- c(min_pen, max_pen)
  n <- length(data)

  test_penalties <- NULL
  numberofchangepoints <- NULL
  penal <- NULL
  overall_cost <- array()
  segmentations <- NULL
  segmentations.smt <- NULL
  #thetas=NULL
  b_between <- array()

  count <- 0

  while (length(pen_interval) > 0){

    new_numcpts <- array()
    new_penalty <- array()
    new_cpts <- array()
    new_smts <- list()
    #new.theta<-array()

    for (b in 1:length(pen_interval)) {
      #ans <-  Fpop(data,pen_interval[b])
      ans <- Rob_seg.std(data, loss="Outlier", lambda=pen_interval[b], lthreshold = lthreshold)
      ## >> GR ADD FOR ROBSEG : compute unpenalized error
      ans$J.est <- ans$cost[length(data)] - pen_interval[b]*length(ans$t.est)
      ## << GR ADD FOR ROBSEG
      resultingcpts <- ans$t.est ##ERROR CORRECTED HERE

      new_numcpts[b] <- length(resultingcpts)-1
      cost.test <- array()
      new_cpts[b]  <- list(resultingcpts)
      new_smts[[b]] <- ans$smt
      # new.theta[b]=list(ans[[5]])
      new_penalty[b] <- ans$J.est
    }

    if (count == 0){
      print(paste("Maximum number of runs of algorithm = ", new_numcpts[1] - new_numcpts[2] + 2, sep = ""))
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }else{
      count <- count + length(new_numcpts)
      print(paste("Completed runs = ", count, sep = ""))
    }

    ## Add the values calculated to the already stored values
    test_penalties <- unique((sort(c(test_penalties,pen_interval))))
    new_numcpts <- c(numberofchangepoints,new_numcpts)
    new_penalty <- c(penal,new_penalty)

    new_cpts <- c(segmentations, new_cpts)
    new_smts <- c(segmentations.smt, new_smts)
    #new.theta <- c(thetas,new.theta)
    numberofchangepoints <- -sort(-new_numcpts) ##can use sort to re-order
    penal <- sort(new_penalty)

    ls <- array()

    for (l in 1:length(new_cpts)){
      ls[l] <- length(new_cpts[[l]])
    }


    ls1 <- sort(ls,index.return = T, decreasing = T)
    ls1 <- ls1$ix


    segmentations <- new_cpts[c(ls1)]
    segmentations.smt <- new_smts[c(ls1)]
    #thetas=new.theta[c(ls1)]

    ## compute new values
    pen_interval <- NULL
    tmppen_interval <- NULL

    for (i in 1:(length(test_penalties)-1)){
      if(abs(numberofchangepoints[i]-numberofchangepoints[i+1])>1){ ##only need to add a beta if difference in cpts>1
        j <- i+1
        tmppen_interval <- (penal[j] - penal[i]) * ((numberofchangepoints[i] - numberofchangepoints[j])^-1)
        pen_interval <- c(pen_interval, tmppen_interval )
      }
    }

    ## discard penalties close to tested one
    if(length(pen_interval)>0){
      for(k in length(pen_interval):1){
        if(min(abs(pen_interval[k]-test_penalties)) < 1e-2) {
          pen_interval=pen_interval[-k]
        }
      }
    }
  }

  ##PRUNE VALUES WITH SAME num_cp
  for(j in length(test_penalties):2){
    if(numberofchangepoints[j]==numberofchangepoints[j-1]){
      numberofchangepoints=numberofchangepoints[-j]
      test_penalties=test_penalties[-j]
      penal=penal[-j]
      segmentations = segmentations[-j]
      segmentations.smt = segmentations.smt[-j]
      #thetas=thetas[-j]
    }
  }



  ###calculate beta intervals
  nb=length(test_penalties)
  beta.int=rep(0,nb)
  beta.e=rep(0,nb)
  for(k in 1:nb){
    if(k==1){
      beta.int[1]=test_penalties[1]
    }else{
      beta.int[k]=beta.e[k-1]
    }
    if(k==nb){
      beta.e[k]=test_penalties[k]
    }else{
      beta.e[k]=(penal[k]-penal[k+1])/(numberofchangepoints[k+1]-numberofchangepoints[k])
    }

  }

  return(list(rbind(test_penalties,beta.int,numberofchangepoints,penal),segmentations, segmentations.smt))
}

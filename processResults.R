directory <- "/n/airoldifs2/lab/afranks/shared_subspace"
lossList <- allLosses <- list()
for(i in 1:3) {
  files <- list.files(directory, pattern=sprintf("srtest-%i",i))
  count <- 0
  lossVecMean <- rep(0, 3)
  lossMat <- matrix(NA, nrow=0, ncol=3)
  for(f in files) {
    print(f)
    count <- count+1
    load(sprintf("%s/%s", directory, f))
    print(resList[[1]]$S)
    lossVecMean <- lossVecMean + lossVec
    lossMat <- rbind(lossMat, lossVec)
  }
  lossList[[i]] <- lossVecMean / count
  allLosses[[i]] <- lossMat
}
apply(allLosses[[3]], 2, function(x) quantile(x, c(0.1, 0.9)))

directory <- "/n/airoldifs2/lab/afranks/shared_subspace"
files <- list.files(directory, pattern="sdim")
lossLst <- c()
for(i in 5:50) {
  print(i)
  load(sprintf("%s/sdimtest-3-%i.RData", directory, i))
  lossLst <- c(lossLst, lossVec[1])
}

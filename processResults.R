directory <- "/n/airoldifs2/lab/afranks/shared_subspace"
lossList <- list()
for(i in 1:3) {
  files <- list.files(directory, pattern=sprintf("vmode-%i",i))
  count <- 0
  lossVecMean <- rep(0, 3)
  for(f in files) {
    print(f)
    count <- count+1
    load(sprintf("%s/%s", directory, f))
    print(resList[[1]]$S)
    lossVecMean <- lossVecMean + lossVec
  }
  lossList[[i]] <- lossVecMean / count

}

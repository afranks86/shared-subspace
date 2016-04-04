directory <- "/n/airoldifs2/lab/afranks/shared_subspace"
for(i in 1:2) {
  files <- list.files(directory, pattern=sprintf("rankTest-%i",i))
  count <- 0
  ssMat <- fullMat <- ratioMat <- matrix(NA, nrow=0, ncol=10)
  for(f in files) {
    print(f)
    count <- count+1
    load(sprintf("%s/%s", directory, f))
    ssMat <- rbind(ssMat, resList[[1]]$ssLoss)
    fullMat <- rbind(fullMat, resList[[1]]$fullLoss)
    ratioMat <- rbind(ratioMat, resList[[1]]$evalRatio)
  }
  save(ssMat, fullMat, ratioMat, file=sprintf("rankResults%i.RData", i))
}



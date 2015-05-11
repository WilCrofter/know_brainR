source("R/voxelize.R")

bigRun <- function(ty1,ty2,nphotons,nreps,ranseed){
 
  phantom <- create_phantom(ty1,ty2)
  f <- function(n){
    result <- get_pairchars(ty1,ty2,phantom,nphotons,ranseed+n)
    #percent exiting through single face, percent absorbed
    c(sum(result$V)/6, sum(result$absorbed))*100/nreps
    
  }
  temp <- sapply(1:nreps,f)
  attr(temp,"ty1") <- ty1
  attr(temp,"ty2") <- ty2 
  attr(temp,"seed") <- ranseed
  temp
}

fullRun <- function(ranseed){
  ans <- list()
  for (ty1 in 1:11){
    for (ty2 in 0:11)
      ans[[length(ans)+1]] <- bigRun(ty1,ty2,100,100,ranseed+100*ty1+ty2)
  }
  ans
}

absorptionMeans <- function(temp){
  lg <- c("Background", "CSF", "Gray Matter", "White Matter", "Fat", "Muscle", "Muscle/Skin", "Skull", "Vessels", "Around fat", "Dura mater", "Bone marrow")
  ans <- data.frame(tissue=lg[2:12],id=1:11,percent.absorption=NA)
  for (i in 1:11){
     ans[i,3] <- mean(sapply(temp[i-1+1:12],function(x){rowMeans(x)})[2,])
  }
  ans
}

flowMeans <- function(temp){
  lg <- c("Background", "CSF", "Gray Matter", "White Matter", "Fat", "Muscle", "Muscle/Skin", "Skull", "Vessels", "Around fat", "Dura mater", "Bone marrow")
  ans <- as.data.frame(matrix(0, 11, 12))
  names(ans) <- lg
  ans <- cbind(tissue=lg[-1], id=1:11, ans)
  for(i in 1:11){
    for(j in 0:11){
      ans[i, j+3] <- mean(temp[[12*(i-1) + j+1]][1,])
    }
  }
  ans
}
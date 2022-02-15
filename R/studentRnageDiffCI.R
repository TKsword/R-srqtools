#' @title Calculate simultaneous confidence intervals
#'
#' @description To obtain simultaneous confidence intervals comparing several binomial parameters using the studentized range distribution with a score statistic.
#'
#' @param name The name of each experiment
#' @param success The number of successful results in each experiment
#' @param size The number of times each experiment
#' @param level The level of Simultaneous confidence intervals
#'
#'
#' @export studentRangeDiffCI
#'
#' @import stats
#'
#' @references Agresti A, Bini M, Bertaccini B, Ryu E. Simultaneous confidence intervals for comparing binomial parameters. Biometrics. 2008;64(4):1270-1275. doi:10.1111/j.1541-0420.2008.00990.x
#'
#' @examples

# 计算多重置信区间

studentRangeDiffCI <- function(name, success, size, level = 0.95){

  no.groups<-length(size)
  no.comp<-choose(no.groups,2)
  Diff.Score.CI.Set<-index<-matrix(c(rep(0,no.comp*3)), ncol=3)
  i<-0
  j<-0
  r<-0
  namelist<-c()
  for(i in 1:(no.groups-1)){
    for(j in (i+1):no.groups){
      r<-r+1
      y1<-success[i]
      y2<-success[j]
      N1<-size[i]
      N2<-size[j]
      conflev<-level
      alpha<-1-level
      diff.int<-diff.student.CI(y1, N1, y2, N2, alpha,no.groups)
      Diff.Score.CI.Set[r,1]<-(y1/N1 - y2/N2)
      Diff.Score.CI.Set[r,2]<-diff.int[1]
      Diff.Score.CI.Set[r,3]<-diff.int[2]
      namelist[r]<-paste0(name[i],"vs",name[j])
    }
  }

  dimnames(Diff.Score.CI.Set)<-list(namelist, c("Estimate", "Lower bound", "Upper bound"))

  return(Diff.Score.CI.Set)
}

#' @noRd
Score.diff<- function (y1,N1,y2,N2,dif){
  p1x<-y1/N1
  p1y<-y2/N2
  nx<-N1
  ny<-N2
  diff = p1x-p1y-dif
  if ( abs(diff) == 0 ) {
    fmdiff = 0
  }
  else{
    t = ny/nx
    a = 1+t
    b = -(1+ t + p1x + t*p1y + dif*(t+2))
    c = dif*dif + dif*(2*p1x + t +1) + p1x + t*p1y
    d = -p1x*dif*(1+dif)
    v = (b/a/3)^3 - b*c/(6*a*a) + d/a/2
    s = sqrt( (b/a/3)^2 - c/a/3)
    if(v>0){
      u=s
    }
    else{
      u=-s
    }
    arg.acos=v/u^3
    if(arg.acos>1){
      arg.acos=1
    }
    w = (3.141592654+acos(arg.acos))/3
    p1d = 2*u*cos(w) - b/a/3
    p2d = p1d - dif
    var = p1d*(1-p1d)/nx + p2d*(1-p2d)/ny
    fmdiff = diff^2/var
  }
  return(fmdiff)
}

#' @noRd
diff.student.CI <- function(y1, N1, y2, N2, alpha, t){
  Q.T.alpha<-qtukey(p=alpha, nmeans=t, lower.tail=FALSE, log.p=FALSE, df=Inf)
  card<-2000
  store.Score.Diff<-c(rep(2,card))
  dif.set<-seq(-.999, .999, length=card )
  k<-1
  for(k in 1:card){
    dif<-dif.set[k]
    score.min<-Score.diff(y1,N1,y2, N2, dif)
    if(score.min < Q.T.alpha^2/2) {store.Score.Diff[k]<-dif}
  }
  left.Score<-min(store.Score.Diff[abs(store.Score.Diff)<2])
  right.Score<-max(store.Score.Diff[abs(store.Score.Diff)<2])
  Score.CI<-c(1,1)*c(left.Score, right.Score)
  return(Score.CI)
}

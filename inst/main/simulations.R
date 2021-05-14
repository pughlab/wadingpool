# libraries 
library(WadingPool)

################################
#### Multiclass Performance ####
N <- 2646

gaugePerf <- function(base_cnt, frac_dists=seq(0.1, 1, by=0.05), seed=1, nstart=10, multi=FALSE){
  if(is.null(names(frac_dists))) names(frac_dists) <- frac_dists
  perfs <- lapply(frac_dists, function(frac_dist){
    cnts <- base_cnt
    cnts[-length(base_cnt)] <- ceiling(cnts[-length(base_cnt)]*frac_dist[-length(base_cnt)])
    
    set.seed(seed)
    draws <- genData(seed = sample(1:10000,1), N=N, cnts = cnts, switch_prob = 0.05)
    hmm1 <- fitHMM(sample='obs', het_cnt = draws, states=length(base_cnt), 
                   nstart = nstart, family=gaussian(), is.sim = TRUE, 
                   ret.raw=TRUE, multi=multi)
    # summary(hmm1$model)
    # boxplot(lapply(split(hmm1$df, hmm1$df$act.state.label), function(i) i$obs))
    # print(c(logLik(hmm1$model), BIC(hmm1$model)))
    
    if(is.null(hmm1)){
      setNames(rep(NA, length(cnts)+1), c(names(cnts), "F1"))
    } else {
      perf <- metrics(hmm1$df$state.label, hmm1$df$act.state.label)
      cbind(t(cnts), "F1"=mean(perf$micro$F1,na.rm=T))
    }
  })
  
  gaus_perfs <- do.call(rbind, perfs)
  rownames(gaus_perfs) <- names(frac_dists)
  return(gaus_perfs)
}

## 2-state Perf
bcnt <- c(400, 50, 10, 5)
cnts_list <- list("large_cnts"=setNames(rep(bcnt[1], 2), c('LOH', 'H1')),
                  "norm_cnts"=setNames(rep(bcnt[2], 2), c('LOH', 'H1')),
                  "low_cnts"=setNames(rep(bcnt[3], 2), c('LOH', 'H1')),
                  "vlow_cnts"=setNames(rep(bcnt[4], 2), c('LOH', 'H1')))
cnts_perf <- lapply(cnts_list, gaugePerf)

cols <- brewer.pal(n = length(cnts_perf)+1, name = "PuBu")
pdf("~/hmm_f1_2state.pdf")
plot(0, type='n', xlim=c(0,1), ylim=c(0,1), xlab="fDist", ylab="F1", las=1)
for(i in seq_along(cnts_perf)){
  perf <- as.data.frame(cnts_perf[[i]])
  points(x=rownames(perf), y=perf$F1, type='l', col=cols[i+1])
  points(x=rownames(perf), y=perf$F1, pch=16, col=cols[i+1])
  axis(side = 1, at = seq(0,1,by=0.2), 
       labels=paste0(round(seq(0,1,by=0.2)*bcnt[i]),
                     "/", bcnt[i]), 
       line = (i*0.5), tick = FALSE, col.axis = cols[i+1], cex.axis=0.7)
}
legend("bottomleft", fill = cols[c(1:length(cnts_perf)+1)], 
       legend=paste0(names(cnts_list), " (", bcnt, ")"), box.lwd = 0)
dev.off()


## 3-state Perf
fdists <- seq(0.1, 1, by=0.1)
fdists <- lapply(setNames(fdists,fdists), function(i) c(i, i+(abs(i-1)/2), 1))
bcnt <- c(400, 50, 10, 5)
cnts_list <- list("large_cnts"=setNames(rep(bcnt[1], 3), c('LOH', 'H1', 'H2')),
                  "norm_cnts"=setNames(rep(bcnt[2], 3), c('LOH', 'H1', 'H2')),
                  "low_cnts"=setNames(rep(bcnt[3], 3), c('LOH', 'H1', 'H2')),
                  "vlow_cnts"=setNames(rep(bcnt[4], 3), c('LOH', 'H1', 'H2')))
cnts_perf <- lapply(cnts_list, gaugePerf, frac_dists=fdists, nstart=20, multi=FALSE)

cols <- brewer.pal(n = length(cnts_perf)+1, name = "PuBu")
pdf("~/hmm_f1_3state.pdf")
plot(0, type='n', xlim=c(0,1), ylim=c(0,1), xlab="fDist", ylab="F1", las=1)
rng <- seq(0,1,by=0.2)
for(i in seq_along(cnts_perf)){
  perf <- as.data.frame(cnts_perf[[i]])
  points(x=as.numeric(rownames(perf)), y=perf$F1, type='l', col=cols[i+1])
  points(x=as.numeric(rownames(perf)), y=perf$F1, pch=16, col=cols[i+1])
  
  xlab <- rng*bcnt[i]
  axis(side = 1, at = rng, 
       labels=paste0(round(xlab), 
                     "/", round(xlab + (abs(bcnt[i]-xlab)/2)), 
                     "/", bcnt[i]), 
       line = (i*0.5), tick = FALSE, col.axis = cols[i+1], cex.axis=0.7)
}
legend("bottomleft", fill = cols[c(1:length(cnts_perf)+1)], 
       legend=paste0(names(cnts_list), " (", bcnt, ")"), box.lwd = 0)
dev.off()

#################################
#### Generate Simulated HMMs ####
fitAndPlot <- function(dat, states){
  for(state in states){
    print(paste0("state = ", state))
    hmm1 <- WadingPool::fitHMM(sample='obs', het_cnt = dat, states=state, family=gaussian(), 
                               is.sim = TRUE, ret.raw=FALSE, multi=FALSE)
    WadingPool::plotHMM(model=hmm1, est_cols = c('state', 'state.label'), 
                        act_cols = c('act.state', 'act.state.label'), xval = 'bin')
  }
}

genData <- function(seed, N, cnts, switch_prob){
  set.seed(seed)
  draws <- WadingPool::simulate(N, cnts = cnts, switch.prob = switch_prob)
  draws$obs <- WadingPool::sanitizeData(draws$obs)
  return(draws)
}


### 2 distributions to sample from
N <- 2646
draws <- genData(seed = 2, N=N, cnts = setNames(c(15, 40), c('LOH', 'H1')), switch_prob = 0.05)
pdf("~/auc_bench_2Act.pdf", height = 10, width = 12)
fitAndPlot(draws, c(1:3))
dev.off()

### 3 distributions to sample from (misfit)
draws <- genData(seed = 1234, N=N, cnts = setNames(c(7, 15, 27), c('LOH', 'H1', 'H2')), switch_prob = 0.05)
pdf("~/auc_bench_3Act.pdf", height = 10, width = 12)
fitAndPlot(draws, c(2:4))
dev.off()


### 5 distributions to sample from (misfit)
draws <- genData(seed = 5, N=N, cnts = setNames(c(5, 15, 30, 40, 45), c('LOH', 'H1', 'H2', 'H3', 'H4')), 
                 switch_prob = 0.05)
pdf("~/auc_bench_5Act.pdf", height = 10, width = 12)
fitAndPlot(draws, c(3:6))
dev.off()


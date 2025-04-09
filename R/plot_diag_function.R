###############
# function to plot Goodness of fit 
# using posterior predictive checks
###############
plot.diag <- function(out, ratio=FALSE, lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  mx <- max(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
  mn <- min(c(samps$dmape.rep[1,], samps$dmape.obs[1,]))
  plot(jitter(samps$dmape.obs[1,]), 
       jitter(samps$dmape.rep[1,]),
       main=paste0("Mean absolute percentage error\nmodel\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$dmape.rep[1,] > samps$dmape.obs[1,]),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot variance/mean ratio
    hist(samps$tvm.rep[1,], nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=samps$tvm.obs[1,1], col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}
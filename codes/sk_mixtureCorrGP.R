rm(list=ls()) 

library('miceadds')
library('gptk')
library('Rfast')
library('graphics')
source.all("./GPtools/util")
source.all("./GPtools")
source('sk_spMultiFullBlockCompute.R')
source('sk_mixtureCorrGPCreate.R')
source('sk_spMultigpComputeM.R')
source('sk_spMultigpParamTrans.R')
source('sk_mixtureCorrGPComputeMat.R')
source('sk_mixtureCorrGPComputeVarDistMat.R')
source('sk_mixtureCorrGPCorLB.R')
source('sk_mixtureCorrGPEstep.R')
source('sk_mixtureCorrGPMstep.R')
source('sk_mixtureCorrGPObjective.R')
source('sk_mixtureCorrGPGradient.R')
source('sk_mixtureCorrGPUpdateParam.R')
source('sk_mixtureCorrGPGradient.R')
source('sk_mixtureCorrGPLLGradient.R')
source('sk_mixtureCorrGPOptimize.R')
source('sk_mixtureCorrGPPriors.R')


###################### Generate data ######################
sf=0.3; or=0.2
obsRangeInput = list(); 
obsRangeInput[[1]] = c(-1,1); obsRangeInput[[2]] = c(-1, 1)
bdata = genData('ggToyMissingPartialLab', diffX=TRUE, sparseFactor=sf, 
                seedVal=200, obsRange=obsRangeInput, obsRate=c(or,or))
X = bdata$X; Y = bdata$y
nout = length(bdata$trueVal$X); bias=bdata$trueVal$bias

colors = c("#FBAD07", "#AB0707")
## Plot data
ymin=min(unlist(bdata$trueVal$ytrue))-2; ymax=max(unlist(bdata$trueVal$ytrue))+2
plot(bdata$X, bdata$y, xlim=c(-1,1), ylim=c(ymin,ymax), 
     type='p', pch=20, col="grey50" );par(new=T)
plot(bdata$X[bdata$observedLab], bdata$y[bdata$observedLab],
     xlim=c(-1,1), ylim=c(ymin,ymax), type='p', 
     pch=20, col=colors[-bdata$lab[bdata$observedLab]+3] );par(new=F)


## Plot data with true trajectories
XTrue = bdata$trueVal$X
yTrue = bdata$trueVal$ytrue
for (m in 1:bdata$nout){
  plot(XTrue[[m]], yTrue[[m]], xlim=c(-1,1), ylim=c(ymin, ymax), 
       type='l', lwd=2, col=colors[-m+3], lty=1); par(new=T)
};par(new=T)
plot(bdata$X, bdata$y, xlim=c(-1,1), ylim=c(ymin,ymax), 
     type='p', pch=20, col=colors[-bdata$lab+3] );par(new=F)

################### Create a MCGP model ###################
priors = mixtureCorrGPPriors(
  N=length(X), nout=nout, labs=bdata$lab, obsLabs=bdata$observedLab
)

nout = length(bdata$trueVal$X)
options = list(nout=nout,
               nlf=1,
               numActive=30,
               beta=1e-5,
               betaTrans='exp',
               fixInducing=TRUE,
               constPriorProb=FALSE,
               constBeta=TRUE,
               priorProb=priors,
               bias=bias,
               obsLab=bdata$observedLab,
               Dir.a = 0.3)

model = mixtureCorrGPCreate(X=X, y=Y, options=options)

####################### Optimization #######################
## initialize parameters
model$param[grep('precision', model$paramName)] = 10
model$param[5] = 100000 

## Optimize the parameters
model = mixtureCorrGPOptimize(model, EM=F, iter=1000, optimizer='SCG')

### Prediction
xPred = bdata$trueVal$X[[1]]

HPrecisionG = model$param[grep('precisionG', model$paramName)]
HSensitivity = model$param[grep('sensitivity', model$paramName)]
HPrecisionU = model$param[grep('precisionU', model$paramName)]
mixtureCorrGPCorLB(model=model, learn='learnall')

yPred = list(); yPredVar = list()
for (m in 1:model$nout){
  Kfstarfstar = ggKernCompute(
    H=c(HPrecisionU, HPrecisionG[m], HSensitivity[m]), xPred, xPred)$K
  Kfstaru = ggXGaussianKernelCompute(
    H = c(HPrecisionG[m],HPrecisionU,HSensitivity[m]), xPred, model$inducingPt)$K
  Dstar = Kfstarfstar - Kfstaru %*% model$Kuuinv %*% t(Kfstaru)
  yPred[[m]] = Kfstaru%*%model$AinvKufDinvy + model$bias[m]
  yPredVar[[m]] = diag(Dstar + Kfstaru%*%model$Ainv%*%t(Kfstaru)) + (1/model$beta[[m]])
}

########################## Plot ##########################
colors = c("#FBAD07", "#AB0707")
colfunc<-colorRampPalette(colors)
par(mar = c(2.1, 2, 1, 1))

ptCol = colfunc(100)[as.numeric(cut(model$varDist$varPi[,1], breaks=100))]

plot(bdata$X, bdata$y, xlim=c(-1,1), ylim=c(ymin, ymax), 
     type='p', pch=20, col='white');par(new=T)
for (m in 1:model$nout){
  polygon(c(xPred, rev(xPred)), 
          c(yPred[[m]]-2.5*yPredVar[[m]], rev(yPred[[m]]+2.5*yPredVar[[m]])), 
          col="grey80", border = NA); par(new=T)
}
for (m in 1:model$nout){
  plot(xPred, yPred[[m]], xlim=c(-1,1), ylim=c(ymin, ymax), 
       type='l', lwd=2, col=colors[-m+3], axes=F); par(new=T)
  plot(XTrue[[m]], yTrue[[m]], xlim=c(-1,1), ylim=c(ymin, ymax), 
       type='l',lty=m+1, axes=F); par(new=T)
}
plot(bdata$X, bdata$y, xlim=c(-1,1), ylim=c(ymin, ymax), type='p', 
     pch=20, col=ptCol, axes=F);par(new=F)

mixtureCorrGPCorLB(model=model, learn='learnall')
mse = rep(NA,model$nout)

for (m in 1:model$nout){ 
  mse[m] = sqrt(sum((bdata$trueVal$y[[m]]-yPred[[m]])^2))/length(yPred[[m]])
}

#probability plot
plot(bdata$X, model$varDist$varPi[,1], pch=20, col=ptCol, ylim=c(-0.1,1.1))

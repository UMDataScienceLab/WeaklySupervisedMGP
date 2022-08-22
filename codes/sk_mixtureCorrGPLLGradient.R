mixtureCorrGPLLGradient <-
  function(model, learn){
    M=model$nout
    
    ### compute partial derivatives with respect to 'covariance' matrix
    # compute dLdKff
    C = model$Ainv + model$AinvKufDinvy %*% t(model$AinvKufDinvy) 
    CKuf = lapply(1:M, function(m){C%*%t(model$Kfu[[m]])})
    KfuAinvKufDinvyy = lapply(1:M, function(m){model$Kfu[[m]]%*%model$AinvKufDinvy%*%model$m[[m]]})
    Dinvy = lapply(1:M, function(m){model$Dinv[[m]]%*%model$m[[m]]})
    Hbar = lapply(1:M, function(m){KfuAinvKufDinvyy[[m]] + t(KfuAinvKufDinvyy[[m]]) - model$Kfu[[m]]%*%CKuf[[m]]})  
    Q = lapply(1:M, function(m){model$Dinv[[m]] - Dinvy[[m]]%*%t(Dinvy[[m]]) + model$Dinv[[m]]%*%Hbar[[m]]%*%model$Dinv[[m]]})
    dLdKff = lapply(1:M, function(m){-0.5*Q[[m]]})  
    
    # compute dLdKuf
    KuuinvKufQ = lapply(1:M, function(m){model$KuuinvKuf[[m]]%*%Q[[m]]})
    CKufDinv = lapply(1:M, function(m){CKuf[[m]]%*%model$Dinv[[m]]})
    AinvKufDinvyyDinv = lapply(1:M, function(m){model$AinvKufDinvy%*%t(Dinvy[[m]])})
    dLdKuf = lapply(1:M, function(m){KuuinvKufQ[[m]] - CKufDinv[[m]] + AinvKufDinvyyDinv[[m]]})
    
    # compute dLdKuu
    QKfuKuuinv = lapply(1:M, function(m){Q[[m]]%*%t(model$KuuinvKuf[[m]])})
    KuuinvKufQKfuKuuinv = Reduce('+',lapply(1:M, function(m){model$KuuinvKuf[[m]]%*%QKfuKuuinv[[m]]}))
    dLdKuu = 0.5*(model$Kuuinv - C - KuuinvKufQKfuKuuinv)
    
    # compute dLdSig
    dLdSig = lapply(1:M, function(m){-0.5*Q[[m]]})
    
    
    ### compute partial derivatives with respect to parameters
    gradParam = rep(0, model$nparam); gradBeta = rep(0, model$nout)
    
    ## compute factors for each parameter and beta
    factorParam = rep(0, model$nparam)
    for (i in 1:model$nparam){
      if (model$nonnegative[i]) factorParam[i] = model$param[i]
      else factorParam[i] = 1
    }
    factorPrecisionG = factorParam[grep("precisionG",model$paramName)]
    factorSensitivity = factorParam[grep("sensitivity",model$paramName)]
    factorPrecisionU = factorParam[grep("precisionU",model$paramName)]
    factorVarLatent = factorParam[grep("variance latent",model$paramName)]
    factorBeta = model$beta 
    
    ################# correction is Needed !! #################
    if (learn == 'learnall'){
      factorVarPiMat = matrix(rep(0,M*(M-1)), ncol=M)
      factorVarPi = list()
      for (n in 1:model$N){
        factorVarPi[[n]] = factorVarPiMat
        for (m1 in 1:M){
          for (m2 in 1:(M-1)){
            if (m1==m2+1) factorVarPi[[n]][m2,m1] = model$varDist$varPi[n,m1]*(1-model$varDist$varPi[n,m1])
            else factorVarPi[[n]][m2,m1] = -model$varDist$varPi[n,(m2+1)]*model$varDist$varPi[n,m1]
          }
        }
      }
    }

    ###########################################################
    
    
    ## compute gradients of parameters
    dLdKffparam = list()
    for (m in 1:M){
      H = c(model$param[[(2*M)+1]], model$param[[m]], model$param[[M+m]])
      dLdKffparam[[m]] = ggKernGradient(H=H, x=model$X[[m]], x2=model$X[[m]], partialMat=dLdKff[[m]])
    }
    
    dLdKufparam = list()
    for (m in 1:M){
      H = c(model$param[[m]], model$param[[(2*M)+1]], model$param[[M+m]])
      dLdKufparam[[m]] = ggXGaussianKernGradient(H=H, x=model$inducingPt, x2=model$X[[m]], partialMat=dLdKuf[[m]])
    }
    
    dLdKuuparam = list()
    for (r in 1:model$nlf){
      H = c(model$param[[(2*M)+r]],1)
      dLdKuuparam[[r]] = rbfKernGradient(H=H, x=model$inducingPt, x2=model$inducingPt, partialMat=dLdKuu)
    }
    
    dLdKuuWhiteparam = list()
    for (r in 1:model$nlf){
      dLdKuuWhiteparam[[r]] = whiteKernGradient(model=model, partialMat=dLdKuu)
    }
    
    # ddLdSigBeta = list()
    # for (m in 1:M){
    #   ddLdSigBeta[[m]] = whiteKernGradient(model=model, partialMat=dLdSig[[m]]) * (-1) * (model$beta[[m]])^(-2)
    # }
    
    ddLdSigBeta = list()
    for (m in 1:M){
      ddLdSigBeta[[m]] = diagKernGradient(model=model, partialMat=dLdSig[[m]], factor=1/as.vector(model$varDist$varPi[,m])) * (-1) * (model$beta[[m]])^(-2)
    }
    
    if (learn == 'learnall'){
      ddLdSigVarPi = list()
      for (m in 1:M){
        ddLdSigVarPi[[m]] = diag(dLdSig[[m]]) * (-1) * (model$varDist$varPi[,m])^(-2) * 1/(model$beta[m])
      }
    }
    
    ######
    # compute dLdprecisionG
    dLdprecisionG = rep(0, M)
    for (m in 1:M){
      for (r in 1:model$nlf){
        dLdprecisionG[m] = dLdprecisionG[m] + dLdKufparam[[m]]$matGradPqr*factorPrecisionG[m]
      }
      dLdprecisionG[m] = dLdprecisionG[m] + dLdKffparam[[m]]$matGradPqr*factorPrecisionG[m]
    }
    
    # compute dLdsensitivity
    dLdsensitivity = rep(0, M)
    for (m in 1:M){
      for (r in 1:model$nlf){
        dLdsensitivity[m] = dLdsensitivity[m] + dLdKufparam[[m]]$gradSenss*factorSensitivity[m]
      }
      dLdsensitivity[m] = dLdsensitivity[m] + dLdKffparam[[m]]$gradSensq*factorSensitivity[m]
    } 
    
    # compute dLdprecisionU
    dLdprecisionU = rep(0, model$nlf)
    for (r in 1: model$nlf){
      dLdprecisionU[r] = dLdprecisionU[r] + dLdKuuparam[[r]]$matGradPrecisionU*factorPrecisionU[r]
      for (m in 1:M){
        dLdprecisionU[r] = dLdprecisionU[r] + dLdKufparam[[m]]$matGradPr*factorPrecisionU[r]
        dLdprecisionU[r] = dLdprecisionU[r] + dLdKffparam[[m]]$matGradPr*factorPrecisionU[r]
      }
    }
    
    # compute dLdlatvar
    dLdlatvar = rep(0, model$nlf)
    for (r in 1:model$nlf){
      dLdlatvar[r] = dLdlatvar[r] + dLdKuuWhiteparam[[r]]*factorVarLatent[r]
    }
    
    # compute dLdbeta
    dLdbeta = rep(0, M)
    for (m in 1:M){
      dLdbeta[m] = dLdbeta[m] + ddLdSigBeta[[m]]*factorBeta[m] + 0.5*sum(model$varDist$varPi[,m]-1)*(1/model$beta[m])*factorBeta[m]
    }
    
    if (model$constBeta) dLdbeta = sum(dLdbeta)
    
    g = c(dLdprecisionG, dLdsensitivity, dLdprecisionU, dLdlatvar, dLdbeta)
    
    if (learn == 'learnall'){ 
      idx.obsT = which(model$obsLab); idx.obsF = which(!model$obsLab)
      
      dLdSigVarPi = matrix(unlist(ddLdSigVarPi), ncol=length(ddLdSigVarPi))
      
      dt2dvarPi = matrix(rep(0, M*model$N), ncol=M)
      dt2dvarPi[idx.obsT,] = -(log(model$varDist$varPi[idx.obsT,]) - log(model$pi[idx.obsT,]) +1)#### 
      digammaPi = digamma(model$Dir.a + model$varDist$varPi) - digamma(M*model$Dir.a + 1)
      dt2dvarPi[idx.obsF,] = -(log(model$varDist$varPi[idx.obsF,]) + 1 - digammaPi[idx.obsF,])
      
      dt3dvarPi = 0.5*(matrix(rep(-log(2*pi*(1/model$beta)),each=model$N),nrow=model$N) - 1/(model$varDist$varPi))
      dLdvarPi = matrix(rep(0, (M-1)*model$N), ncol=(M-1))
      # dKLDirdvarPi = -log()
      for (n in 1:model$N){
        # dLdvarPi[n] = (dLdSigVarPi[n,]+dt3dvarPi[n,]) %*% t(factorVarPi[[n]])
        dLdvarPi[n] = (dLdSigVarPi[n,]+dt2dvarPi[n,] + dt3dvarPi[n,]) %*% t(factorVarPi[[n]])
        # dLdvarPi[n] = (dt3dvarPi[n,]) %*% t(factorVarPi[[n]])
      }
      # dLdSigVarPi
      # for (m in 2:M){  ##### m is indeed m-1 !!
      #   dLdvarPi[,(m-1)] = (dt3dvarPi[,1])*factorVarPi[,1] + (dt3dvarPi[,m])*factorVarPi[,m]
      # }
      # # dLdvarPi[1] = 1/model$varDist$varPi[1,1] * factorVarPi[1,1] 
      dLdvarPi = as.vector(dLdvarPi)
      g = c(dLdprecisionG, dLdsensitivity, dLdprecisionU, dLdlatvar, dLdbeta, dLdvarPi)
    }
    
    
    
    return(g)
  }
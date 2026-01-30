# A function for performing the Instrumental Variables approach by Kelejian and Prucha (1998).



IV <- function(nam, ...){
  IV.Effects.Freq = function(y,X,A,alpha=0.05,alternative='twoSided'){
    ZZ = cbind(X,A%*%y)
    if(mean(X[,1] == 1) == 1){
      X1 = X[,-1]
    }else{
      X1 = X
    }
    AX = A%*%X1
    HH = model.matrix(~X1+AX+A%*%AX)
    # HH = cbind(X,AX,A%*%AX)
    # PP = tcrossprod(HH%*%qr.solve(crossprod(HH)),HH)
    PP = tcrossprod(HH%*%chol2inv(chol(crossprod(HH))),HH)
    ZHat = PP%*%ZZ
    ZtZInv= qr.solve(crossprod(ZHat)) 
    # ZtZInv= qr.solve(crossprod(ZHat,ZZ)) #Mistake in the sphet vignette.  Checked with KP98
    
    Coefs = ZtZInv%*%crossprod(ZHat,y)
    
    s2Hat = drop(crossprod(y-ZZ%*%Coefs))/nrow(X)
    covMat = s2Hat*ZtZInv
    
    Coefs = data.frame(Estimates=Coefs,
                       SE=sqrt(diag(covMat)))
    Coefs = data.frame(Coefs,
                       LB=Coefs$Est+qnorm(alpha/2)*Coefs$SE,
                       UB=Coefs$Est+qnorm(1-alpha/2)*Coefs$SE,
                       pvalue=2^(alternative=='twoSided')*
                         pnorm(-abs(Coefs$Est/Coefs$SE)))
    rownames(Coefs)[nrow(Coefs)] = 'networkEffect'
    return(list(coefs=Coefs,s2=s2Hat,covMat=covMat))
  }
  
  # This is unfinished
  IV.Disturb.Freq = function(y,X,A,alpha=0.05,alternative="twoSided"){
    
  }
  
}

source("MDPDE.r", local = TRUE)
MDPDE_BETA2 <- function(y, X, Z, q=0.6){
  
  truncated <- 1e-9
  Y<-y
  Y[Y <= truncated] = truncated
  Y[Y >= 1-truncated] = 1-truncated
  y <- Y # response value
  
  error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                      linkmu="logit", linkphi="identity", weights=FALSE))
  

  if(class(error)!="try-error"){
    return(fit_MDPDE)
    
  } else {
    truncated1 <- truncated*10
    Y[Y <= truncated1]=truncated1
    Y[Y >= 1 - truncated1]= 1 - truncated1
    y <- Y # response value
    
    error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                        linkmu="logit", linkphi="identity", weights=FALSE))
    if(class(error)!="try-error"){
      return(fit_MDPDE)
      
    } else {
      truncated2 <- truncated1*10
      Y[Y <= truncated2]=truncated2
      Y[Y >= 1 - truncated2]= 1 - truncated2
      y <- Y # response value
      
      
      error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                          linkmu="logit", linkphi="identity", weights=FALSE))
      
      if(class(error)!="try-error"){
        return(fit_MDPDE)
        
      } else {
        truncated3 <- truncated2*10
        Y[Y <= truncated3]=truncated3
        Y[Y >= 1 - truncated3]= 1 - truncated3
        y <- Y # response value
        
        error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                            linkmu="logit", linkphi="identity", weights=FALSE))
        
        if(class(error)!="try-error"){
          return(fit_MDPDE)
          
        } else {
          truncated4 <- truncated3*10
          Y[Y <= truncated4]=truncated4
          Y[Y >= 1 - truncated4]= 1 - truncated4
          y <- Y # response value
          
          error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                              linkmu="logit", linkphi="identity", weights=FALSE))
          
          if(class(error)!="try-error"){
            return(fit_MDPDE)
            
          } else {
            truncated5 <- truncated4*10
            Y[Y <= truncated5]=truncated5
            Y[Y >= 1 - truncated5]= 1 - truncated5
            y <- Y # response value
            
            error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                                linkmu="logit", linkphi="identity", weights=FALSE))
            
            if(class(error)!="try-error"){
              return(fit_MDPDE)
              
            } else {
              truncated6 <- truncated5*10
              Y[Y <= truncated6]=truncated6
              Y[Y >= 1 - truncated6]= 1 - truncated6
              y <- Y # response value
              
              error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                                  linkmu="logit", linkphi="identity", weights=FALSE))
              
              if(class(error)!="try-error"){
                return(fit_MDPDE)
                
              } else {
                Y<-(Yy*(length(Yy[,1])-1)+0.5)/length(Yy[,1])
                y <- Y # response value
                
                error<- try(fit_MDPDE <- MDPDE_BETA(y=y, X=X, Z=Z, qoptimal=TRUE, q0=q, m=3, L=0.02, qmin=0.5, spac=0.02, method="BFGS", startV="CP",
                                                    linkmu="logit", linkphi="identity", weights=FALSE))
                
                if(class(error)!="try-error"){
                  return(fit_MDPDE)
                  
                } else {
                  next
                }
              }
            }
          } 
        }
      }
    }
  }
  
  
  
}
  
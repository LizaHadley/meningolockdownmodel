
meningo_fun_novax <- function(t,x,params){  # Pre-vaccination equations (pre-2015)
  nage <- 100
  
  sf <- 1
  nf <- 2 
  mf <- 3
  inc <- 4
  ninc <- 5
  
  S0 <- x[index2(sf,1,nage)]  
  S <- x[index2(sf,2:nage,nage)]
  N0 <- x[index2(nf,1,nage)]
  N <- x[index2(nf,2:nage,nage)]
  M0 <- x[index2(mf,1,nage)]
  M <- x[index2(mf,2:nage,nage)]
  Inc0 <- x[index2(inc,1,nage)]
  Inc <- x[index2(inc,2:nage,nage)]
  NInc0 <- x[index2(ninc,1,nage)]
  NInc <- x[index2(ninc,2:nage,nage)]
  
  print(t)
  
  with(as.list(params),{
    
    for(i in 1:100){
      lambda_n[i] <- beta[i]*sum(L[i,]*prev)  # adding social mixing in N compartments (to allow exploration of pandemic effects on all strains)
    }
    
    lambda_n <- lambda_n*(1-p.vac)
    
    dS0 <-   -lambda_n[1]*S0  - beta[1]*(sum(as.matrix(L[1,-1])*as.matrix(M)))*S0 - beta[1]*(as.matrix(L[1,1]*as.matrix(M0)))*S0 + r*(N0+M0) - mort[1]*S0 + (ageing[1]+mort[1])*P[1] - ageing[1]*S0
    
    dS <- -lambda_n[-1]*S - beta[-1]*(as.matrix(L[-1,-1])%*%as.matrix(M))*S - beta[-1]*(as.matrix(L[-1,1])%*%as.matrix(M0))*S + r*(M+N) - mort[-1]*S + ageing[-nage]*c(S0,S[-nage+1]) - ageing[-1]*S  
    
    dN0 <- lambda_n[1]*S0 - r*N0 - mort[1]*N0 - ageing[1]*N0
    
    dN <- lambda_n[-1]*S - r*N -  mort[-1]*N + ageing[-nage]*c(N0,N[-(nage-1)]) - ageing[-1]*N               
    
    dM0 <-  - r*M0 + beta[1]*(sum(as.matrix(L[1,-1])*as.matrix(M)))*S0 + beta[1]*(as.matrix(L[1,1]*as.matrix(M0)))*S0 - mort[1]*M0 - ageing[1]*M0 
    
    dM <- -r*M + beta[-1]*(as.matrix(L[-1,-1])%*%as.matrix(M))*S + beta[-1]*(as.matrix(L[-1,1])%*%as.matrix(M0))*S - mort[-1]*M + ageing[-nage]*c(c(M0,M[-(nage-1)])) - ageing[-1]*M
    
    dInc0 <- beta[1]*(sum(as.matrix(L[1,-1])*as.matrix(M)))*S0 + beta[1]*(as.matrix(L[1,1]*as.matrix(M0)))*S0
    
    dInc <- beta[-1]*(as.matrix(L[-1,-1])%*%as.matrix(M))*S + beta[-1]*(as.matrix(L[-1,1])%*%as.matrix(M0))*S  # tracking incidence: same as dM but only transmission terms (no death, ageing etc)
    
    dNInc0 <- lambda_n[1]*S0
    
    dNInc <- lambda_n[-1]*S
    
    out <- c(dS0,dS,dN0,dN,dM0,dM,dInc0,dInc,dNInc0,dNInc)
    list(out) 
  })
}







meningo_fun_vaccination <- function(t,x,params){  # Vaccination equations (from 2015)
  nage <- 100
  
  sf <- 1
  nf <- 2 
  mf <- 3
  vsif <- 4
  vnif <- 5
  vmif <- 6
  inc <- 7
  vinc <- 8
  ninc <- 9
  
  S0 <- x[index2(sf,1,nage)]  
  S <- x[index2(sf,2:nage,nage)]
  N0 <- x[index2(nf,1,nage)]
  N <- x[index2(nf,2:nage,nage)]
  M0 <- x[index2(mf,1,nage)]
  M <- x[index2(mf,2:nage,nage)]
  VSI0 <- x[index2(vsif,1,nage)]
  VSI <- x[index2(vsif,2:nage,nage)]
  VNI0 <- x[index2(vnif,1,nage)]
  VNI <- x[index2(vnif,2:nage,nage)]
  VMI0 <- x[index2(vmif,1,nage)]
  VMI <- x[index2(vmif,2:nage,nage)]
  Inc0 <- x[index2(inc,1,nage)]
  Inc <- x[index2(inc,2:nage,nage)]
  VInc0 <- x[index2(vinc,1,nage)]
  VInc <- x[index2(vinc,2:nage,nage)]
  NInc0 <- x[index2(ninc,1,nage)]
  NInc <- x[index2(ninc,2:nage,nage)]
  
  
  print(t)
  with(as.list(params),{
    
    Lt <- (1-rd_contact(t/365,rd_drop1,rd_drop2,rd_duration))*L
    
    for(i in 1:100){
      lambda_n[i] <- beta[i]*sum(Lt[i,]*prev)
    }
    
    lambda_n <- lambda_n*(1-p.vac)
    
    dS0 <-   -lambda_n[1]*S0  - beta[1]*(sum(as.matrix(Lt[1,-1])*as.matrix(M+VMI)))*S0 - beta[1]*(as.matrix(Lt[1,1]*as.matrix(M0+VMI0)))*S0 + r*(N0+M0) - mort[1]*S0 + (ageing[1]+mort[1])*P[1] - ageing[1]*S0
    
    dS <- -lambda_n[-1]*S - beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*S - beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*S + r*(M+N) - mort[-1]*S + ageing[-nage]*c(S0,S[-nage+1]) - ageing[-1]*S  + w*VSI - rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(S0,S[-nage+1])
    
    dN0 <- lambda_n[1]*S0 - r*N0 - mort[1]*N0 - ageing[1]*N0
    
    dN <- lambda_n[-1]*S - r*N -  mort[-1]*N + ageing[-nage]*c(N0,N[-(nage-1)]) - ageing[-1]*N + w*VNI  - rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(N0,N[-nage+1])      
    
    dM0 <-  - r*M0 + beta[1]*(sum(as.matrix(Lt[1,-1])*as.matrix(M+VMI)))*S0 + beta[1]*(as.matrix(Lt[1,1]*as.matrix(M0+VMI0)))*S0 - mort[1]*M0 - ageing[1]*M0 
    
    dM <- -r*M + beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*S + beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*S - mort[-1]*M + ageing[-nage]*c(M0,M[-(nage-1)]) - ageing[-1]*M + w*VMI - rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(M0,M[-nage+1])
    
    dVSI0 <-  0  # no babies receive the MenACWY vaccine
    
    dVSI <- rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(S0,S[-nage+1]) + ageing[-nage]*c(VSI0,VSI[-nage+1]) - ageing[-1]*VSI  - mort[-1]*VSI + r*(VNI+VMI) -  lambda_n[-1]*VSI - (1-kappa)*beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*VSI - (1-kappa)*beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*VSI - w*VSI
    
    dVNI0 <- 0
    
    dVNI <- rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(N0,N[-nage+1]) + ageing[-nage]*c(VNI0,VNI[-nage+1]) - ageing[-1]*VNI - mort[-1]*VNI - r*VNI + lambda_n[-1]*VSI - w*VNI
    
    dVMI0 <- 0
    
    dVMI <- rate_1y(t/365,drop)*rate_vax*ageing[-nage]*c(M0,M[-nage+1]) + ageing[-nage]*c(VMI0,VMI[-nage+1]) - ageing[-1]*VMI - mort[-1]*VMI - r*VMI + (1-kappa)*beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*VSI  + (1-kappa)*beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*VSI - w*VMI
    
    dInc0 <- beta[1]*(sum(as.matrix(Lt[1,-1])*as.matrix(M+VMI)))*S0 + beta[1]*(as.matrix(Lt[1,1]*as.matrix(M0+VMI0)))*S0
    
    dInc <- beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*S + beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*S  # tracking incidence: same as dM but only transmission terms (no death, ageing etc)
    
    dVInc0 <- 0
    
    dVInc <- (1-kappa)*beta[-1]*(as.matrix(Lt[-1,-1])%*%as.matrix(M+VMI))*VSI  + (1-kappa)*beta[-1]*(as.matrix(Lt[-1,1])%*%as.matrix(M0+VMI0))*VSI
    
    dNInc0 <- lambda_n[1]*S0 
    
    dNInc <- lambda_n[-1]*S + lambda_n[-1]*VSI
    
    out <- c(dS0,dS,dN0,dN,dM0,dM,dVSI0,dVSI,dVNI0,dVNI,dVMI0,dVMI,dInc0,dInc,dVInc0,dVInc,dNInc0,dNInc)
    list(out,Lt[1,1]) 
  })
}

# Inc = incidence (NEW carriage) for M compartment
# VInc = incidence (NEW carriage) for VMI compartment
# NInc = incidence (NEW carriage) for N and VNI compartments
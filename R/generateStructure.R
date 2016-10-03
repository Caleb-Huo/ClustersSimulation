generateStructure <- function(anModule,Gmean,Gsd,aSindex,sigma1,sigma2,rho = 0.5,df.prior = 60,meanGenesPerModule=20,
					amuLayer1=NULL,anGClust=NULL,constrain=FALSE,diffmu=NULL,fold=1){
  aSubtypes = length(aSindex)
  ## generate uniform template for each template
  ## mean 9,sd 2
  if(is.null(amuLayer1)){
      amuLayer1 = matrix(rnorm(anModule*aSubtypes,Gmean,Gsd),nrow=anModule,ncol=aSubtypes)
  	
	  while(constrain){
	  	if(is.null(diffmu)){
	  		stop("please specify diffmu, which is a constrain parameter!")
	  	}
	  	if(sd(amuLayer1)<diffmu){
	  		stop("sd(amuLayer1)<diffmu, it is unlikely to statisfy this constrain")
	  	}
	  	resampleIndex = apply(amuLayer1,1,function(x) max(x)-min(x))<diffmu
	  	resampleLen = sum(resampleIndex)
	  	if(resampleLen==0){
	  		constrain=FALSE
	  		next
	  	}
	  	amuLayer1[resampleIndex,] = matrix(rnorm(resampleLen*aSubtypes,Gmean,Gsd),nrow=resampleLen,ncol=aSubtypes)
	  }
  }

  if(is.null(anGClust)){
    anGClust = rep(meanGenesPerModule, anModule)        
  }

  an = sum(sapply(aSindex,length))
  amuLayer2  = matrix(NA,nrow=nrow(amuLayer1),ncol=an)  
  ## expand the gene template for the sample direction
  
  amuLayer1_fold <- t(apply(amuLayer1,1,function(x) (x-min(x))*fold + min(x) ))
  print(amuLayer1_fold)
  for(i in 1:nrow(amuLayer1)){
    for(j in 1:aSubtypes){
      tmpIndex = aSindex[[j]]
      amuLayer2[i,tmpIndex] = rnorm(n=length(tmpIndex),mean=amuLayer1_fold[i,j],sd=sigma1)
    }
  }
  dim(amuLayer2)
  #gplots::heatmap.2(amuLayer2, col="greenred" ,trace="none",Rowv=NA,  Colv=NA,keysize=1.3)  

  ## add correlated gene structure
  aresData1 = matrix(NA,nrow = sum(anGClust), ncol = an)
  subtypeModuleIndex <- NULL
  for(i in 1:anModule){
    if(i==1){
      firstIndex=1
    }else{
      firstIndex = cumsum(anGClust)[i-1]+1    
    }    
    endIndex = cumsum(anGClust)[i]
    tmpRowIndex = firstIndex:endIndex
	subtypeModuleIndex[[i]] <- tmpRowIndex
    aresData1[tmpRowIndex,] = matrix(rep(amuLayer2[i,],anGClust[i]),nrow=anGClust[i],byrow=TRUE)+sample.correlated.data(anGClust[i], an, rho, df.prior, sigma2)
  }
  
  #gplots::heatmap.2(aresData1, col="greenred" ,trace="none",Rowv=TRUE,  Colv=NA,keysize=1.3)      
  #range(aresData1)
  
  res = list(data=aresData1,subtypeModuleIndex=subtypeModuleIndex)
  return(res)
}

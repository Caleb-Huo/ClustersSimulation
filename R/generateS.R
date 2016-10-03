##' generate S for ISKmeans
##'
##' generate S for ISKmeans
##' @title generateS
##' @param seed random seed
##' @param S number of studies
##' @param k number of clusters
##' @param meanSamplesPerK mean samples per cluster
##' @param nModule number of modules. A module is a group of genes.
##' @param meanGenesPerModule number of genes per module
##' @param Gmean gene expression template follows N(Gmean,Gsd^2)
##' @param Gsd gene expression template follows N(Gmean,Gsd^2)
##' @param sigma1 noise 1
##' @param sigma2 noise 2
##' @param sigma3 noise 3
##' @param G0 number of noise genes
##' @param nconfounder number of confounders
##' @param nrModule number of modules for confounding variables
##' @param rMeanSubtypes number of subtypes defined by confounding variables
##' @param diffmu effect size difference for subtype predictive genes
##' @param fold how to vary subtype predictive gene signal. 1: original. 0: no signal.
##' @param rho para for inverse Wishart distribution.
##' @param df.prior para for inverse Wishart distribution.
##' @param groupProb subtype predictive genes have prior group information. By prob 1-groupProb, the information will be altered.
##' @return alist
##' @author Caleb
##' @export
##' @examples
##' Sdata =generateS(seed=15213,S=2,k=3,meanSamplesPerK=c(40,40,30),nModule=30,meanGenesPerModule=30,
##'                  sigma1=1,sigma2=1,sigma3=1,G0=5000,
##'                   nconfounder=4,nrModule=20,rMeanSubtypes=3,diffmu=1,fold=c(1,1),
##' 				  rho = 0.5,df.prior = 100,groupProb=1)
##' #
##' Sdata$subPredictGeneUnion
##' sum(Sdata$subPredictGeneUnion)
generateS <- function(seed=15213,S=2,k=3,meanSamplesPerK,nModule,meanGenesPerModule=20,
                      Gmean=9,Gsd=2,sigma1,sigma2,sigma3,G0,
                      nconfounder,nrModule,rMeanSubtypes,diffmu,fold=rep(1,S),
                      rho = 0.5,df.prior = 60,groupProb=1){
  ## seed: random number
  ## S: number of studies
  ## k: number of main subtypes
  ## meanSamplesPerK: mean number of Samples Per subtype (K)
  ## nModule: number of gene clusters (gene Module)
  ## meanGenesPerModule: mean number of Genes Per gene Module
  ## mmin: min parameter from the UNIF distributin
  ## mmax: max parameter from the UNIF distributin
  ## cmax: max parameter from the UNIF distributin, for the confounding variable.
  ## sigma1: standard deviation of the within subtype (including confounders)
  ## sigma2: standard deviation of inverse Wischart covariate matrix
  ## sigma3: standard deviation of the house keeping genes
  ## G0: number of house keeping genes
  ## nconfounder: number of confounders
  ## nrModule: mean number of gene modules for each confounder
  ## rMeanSubtypes: number of confoundering subtypes
  ## rho: parameter used in the inverse Wischart structure
  ## df.prior: parameter used in the inverse Wischart structure

  ## when compare with TSS and scale SKM,

  ## within cluster variance first, then correlated genes
  set.seed(seed)

  ## main function:
  ## generate S study
  result=NULL
  ## gene template for the main part
  muLayer1=NULL
  ## number of correlated genes for the main part
  nGClust=NULL

  ## gene template for the confounding part
  ## number of correlated genes for the confounding part
  rnGClust = NULL

  ## determine number of subtypes for the confounding variable
  rk <- rep(k,nconfounder)

  ## prepare samples for the subtypes
  nall = meanSamplesPerK
  cumnall = cumsum(nall)
  n = sum(nall)
  Sindex = NULL
  label = numeric(n)
  for(i in 1:length(nall)){
    if(i==1){
      Sindex[[i]] = 1:cumnall[i]
    } else {
      Sindex[[i]] = (cumnall[i-1]+1):cumnall[i]
    }
    ## label the samples
    label[Sindex[[i]]] = i
  }


  for(s in 1:S){
    dataMain = generateStructure(nModule,Gmean,Gsd,Sindex,sigma1,sigma2,rho = rho,
		df.prior = df.prior,meanGenesPerModule=meanGenesPerModule,
		amuLayer1=muLayer1,anGClust=nGClust,constrain=TRUE,diffmu=diffmu,fold=fold[s])


    resData = dataMain$data
    subtypeModuleIndex = dataMain$subtypeModuleIndex
    #gplots::heatmap.2(resData, col="greenred" ,trace="none",Rowv=TRUE,  Colv=NA,keysize=1.3)

    ## next step, simulate correlated random genes, with rk~POI(k) clusters
    ## the biology plausibility is that there are other factors would influnce the gene
    ## expression, e.g. gender, race, geological information, disease grade, level of hormone
    ## then exp~UNIF(), for each of the rk subclass
    ## generate m unrelated subtypes, for each m, there should be mmoduls

    if(length(rk)!=0){
	    for (i in 1:length(rk)){
	      rSindex = rPartationSample(n,rk[i])
	  	  rdatai = generateStructure(nrModule,Gmean,Gsd,rSindex,sigma1,sigma2,rho = rho,df.prior = df.prior,meanGenesPerModule=meanGenesPerModule)
	      resData = rbind(resData,rdatai$data)
	    }
    }
    dim(resData)
	confounderIndex <- (nrow(dataMain$data) + 1):nrow(resData)

    #gplots::heatmap.2(resData, col="greenred" ,trace="none",Rowv=TRUE,  Colv=NA,keysize=1.3)

    ## finally, add some house keeping genes
    ## sigma3 = 1, error of the random genes
    templateRandom = rnorm(G0,Gmean,Gsd)

    randomPart  = matrix(NA,nrow=G0,ncol=n)  ## 20*n
    for(i in 1:G0){
      randomPart[i,] = rnorm(n=n,mean=templateRandom[i],sd=sigma3)
    }

    resData = rbind(resData,randomPart)

	randomIndex <- (nrow(resData) - nrow(randomPart) + 1):nrow(resData)

    result[[s]]=list(S=resData,label=label,subtypeModuleIndex=subtypeModuleIndex,confounderIndex=confounderIndex,randomIndex=randomIndex)
  }

  group <- NULL

  for(s in 1:length(result)){
	  aresult <- result[[s]]
	  agroup <- aresult$subtypeModuleIndex
	  noiseIndex <- c(aresult$confounderIndex, aresult$randomIndex)
	  for(g in 1:length(agroup)){
		  bgroup <- agroup[[g]]
		  relinkIndex <- rbinom(length(bgroup),1,groupProb) == 0
		  bgroup[relinkIndex] = sample(noiseIndex, sum(relinkIndex), replace=FALSE)
		  agroup[[g]] <- bgroup
	  }
	  group[[s]] <- agroup
  }

  omics <- lapply(result,function(x) x$S)
  labels <- result[[1]]$label
  subPredictGene <- lapply(result,function(x) unlist(x$subtypeModuleIndex))
  subPredictGeneUnion <- NULL
  for(s in 1:length(result)){
	  ap <- nrow(omics[[s]])
	if(s==1){
		subPredictGeneUnion <- 1:ap %in% subPredictGene[[s]]
		groupUnion <- group[[s]]
		d <- t(omics[[s]])
	} else {
	  	subPredictGeneUnion <- c(subPredictGeneUnion, 1:ap %in% subPredictGene[[s]])
		agroup <- group[[s]]
		for(g in 1:length(agroup)){
			groupUnion[[g]] <- c(groupUnion[[g]], (agroup[[g]] + ap))
		}
		d <- cbind(d,t(omics[[s]]))

	}
  }

  finalResult <- list(omics=omics, d=d, labels=labels, group=group, subPredictGene=subPredictGene, groupUnion=groupUnion, subPredictGeneUnion=subPredictGeneUnion)

  return(finalResult)
}

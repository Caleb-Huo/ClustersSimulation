rPartationSample <- function(n,rSubtypes){
  rSindex = vector('list',rSubtypes)
  meanSample = floor(n/rSubtypes)
  samplePool = 1:n
  for(i in 1:rSubtypes){
    if(i==rSubtypes){
      rSindex[[i]] = samplePool
    } else{
	  aSampleNum = rpois(1,meanSample)
	  if(aSampleNum>=length(samplePool) | length(samplePool)<n/rSubtypes/3){
		 return(rPartationSample(n,rSubtypes))
	  }
      tmpGroup = sample(samplePool,aSampleNum)
      rSindex[[i]] = tmpGroup
      samplePool = setdiff(samplePool,tmpGroup)
    }
  }
  return(rSindex)
}
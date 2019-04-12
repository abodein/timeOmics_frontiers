### variables with a non 0 value on the components # (3 types of data in that case)
summary.block <- function(block){
  
  sum_load = NULL
  res = NULL
  
  for (i in 1:3) {
    
    #for each block, loading of the variables on the components is checked. Only variables with loading different from 0 are kept.
    sum_load[[i]]=apply(block$loadings[[i]],1,sum)
    res[[i]] = which(sum_load[[i]]!=0)
  
    }
  return(list(summary.OTU=res[[1]],summary.metabolites=res[[2]],summary.perf=res[[3]]))
}


### components where variables have a non 0 loading
compo.block <- function(res.block){
  
  summary = summary.block(res.block)
  
  test_null = NULL
  no_compo = NULL
  sign_compo = NULL
  
  no_compo_total = NULL
  sign_compo_total = NULL
  
  for (i in 1:3){

    
    #select variables with a non 0 loading 
    test_null[[i]] = res.block$loadings[[i]][summary[[i]],]
    
    #for all the components, all the variables
    for (j in 1:dim(test_null[[i]])[1]) {

      
      for (k in 1:dim(test_null[[i]])[2]) {

        
        #if variable has a non 0 loading on comp k
        if (test_null[[i]][j,k] != 0) {
          
          #number of the component 
          no_compo[j] = k
          
          #sign_compo : positive or negative loading on that component 
          sign_compo[j] = sign(test_null[[i]][j,k])
          
        }
        
        
      }
      
    }
    
    #name and sign of the components for each  block
    no_compo_total[[i]]= no_compo
    sign_compo_total[[i]] = sign_compo
    
    #give a name to the groups
    sign_compo_total[[i]]= str_replace(as.character(sign_compo_total[[i]]), "-1", "neg" )
    sign_compo_total[[i]] = str_replace(sign_compo_total[[i]], "1", "pos" )
    
    names(no_compo_total[[i]])=rownames(test_null[[i]])
    names(sign_compo_total[[i]])=rownames(test_null[[i]])
    
    no_compo = NULL
    sign_compo = NULL
    
  }
  
  return(list( var.selected.OTU = no_compo_total[[1]], var.selected.metabolites = no_compo_total[[2]], 
               var.selected.perf = no_compo_total[[3]],
               sign.OTU = sign_compo_total[[1]], sign.metabolites = sign_compo_total[[2]], 
               sign.perf = sign_compo_total[[3]]))
  
}





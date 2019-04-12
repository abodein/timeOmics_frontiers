# useful functions

# we separate the positively from the negatively correlated profiles then calculate the mean expression values per 'cluster' (comp) on the original mat data
# we save the selected variables, which comp they were selected and their pos / neg sign

spls.result = function(data,     # the data on which to calculate the mean profiles
                       block,    # either 'X' or 'Y'
                       spls.res  # result from a sPLS object
                       ){
  
  # outputs:
  # - the selected variables and on which component
  # - the sign of each variable 
  # - the mean profiles per pos / negative profiles per component

  
  
  # initialise
spls.means = NULL
var.selected = var.sign = rep(NA, ncol(data))
names(var.selected) = names(var.sign) = colnames(data)
profile.old = character()

for(comp in 1:spls.res$ncomp){
  # extract selected variables, depending on the block
  if(block == 'X'){
    profile = selectVar(spls.res, comp = comp)$X$name
    index.pos = which(selectVar(spls.res, comp = comp)$X$value >= 0)
  }else{      # if block is 'Y'
    profile = selectVar(spls.res, comp = comp)$Y$name
    index.pos = which(selectVar(spls.res, comp = comp)$Y$value >= 0)
  }
  # name of feature with positive weight
  profile.pos = profile[index.pos]
  # name of feature with negative weight
  profile.neg = setdiff(profile, profile.pos)
  
  # fill the means matrix only if there are orthogonal features in the selection
  # store where the feature was selected, and whether pos / neg sign
  var.selected[setdiff(profile, profile.old)] = comp
  name.pos = setdiff(profile.pos, profile.old)
  name.neg = setdiff(profile.neg, profile.old)
  
  if(length(name.pos) != 0 & length(name.neg) != 0){
    var.sign[name.pos] = 'pos'
    var.sign[name.neg] = 'neg'
    mean.tmp = t(as.matrix(apply(data[,setdiff(profile, profile.old)],1,function(x)tapply(x, var.sign[setdiff(profile, profile.old)], mean))))
    colnames(mean.tmp) = paste0('comp', comp, colnames(mean.tmp)) 
    spls.means = cbind(spls.means, mean.tmp)
  }else if(length(name.pos) != 0){
    var.sign[name.pos] = 'pos'
    spls.means = cbind(spls.means,  apply(as.matrix(data[, name.pos]), 1, mean))
    colnames(spls.means) = c(colnames(spls.means)[-ncol(spls.means)], paste0('comp', comp, 'pos'))
  }else{  # if length(name.pos) = 0
    var.sign[name.neg] = 'neg'
    spls.means = cbind(spls.means,  apply(as.matrix(data[, name.neg]), 1, mean))
    colnames(spls.means) = c(colnames(spls.means)[-ncol(spls.means)], paste0('comp', comp, 'neg'))
  }
  profile.old = c(profile.old, profile)
}

# this shows the mean profile plots within the pos / neg profiles, per component
# careful the scale can be misleading!
col.remove = unique(which(is.na(spls.means), arr.ind = TRUE)[, 'col'])
if(length(col.remove) >0) spls.means = spls.means[, -c(col.remove)]


# all profile plots, on the original variables, selected
# remove variables not selected
var.selected = var.selected[!is.na(var.selected)]
#length(var.selected)   
var.sign = var.sign[!is.na(var.sign)]
#length(var.sign)

return(list(var.selected = var.selected, var.sign = var.sign, spls.means = spls.means))

}

diablo.result = function(data,     # the data on which to calculate the mean profiles
                         block,    # bloc name
                         diablo.res  # result from a diablo object
){
  
  # outputs:
  # - the selected variables and on which component
  # - the sign of each variable 
  # - the mean profiles per pos / negative profiles per component
  
  
  
  # initialise
  diablo.means = NULL
  var.selected = var.sign = rep(NA, ncol(data))
  names(var.selected) = names(var.sign) = colnames(data)
  profile.old = character()
  
  for(comp in 1:diablo.res$ncomp){
    # eOTUtract selected variables, depending on the block
    if(block == 'OTU'){
    profile = selectVar(spls.res, comp = comp)$OTU$name
    index.pos = which(selectVar(diablo.res, comp = comp)$OTU$value >= 0)
  }else{      # if block is 'metabo'
    profile = selectVar(diablo.res, comp = comp)$metabo$name
    index.pos = which(selectVar(diablo.res, comp = comp)$metabo$value >= 0)
  }
  # name of feature with positive weight
  profile.pos = profile[index.pos]
  # name of feature with negative weight
  profile.neg = setdiff(profile, profile.pos)
  
  # fill the means matriOTU onlmetabo if there are orthogonal features in the selection
  # store where the feature was selected, and whether pos / neg sign
  var.selected[setdiff(profile, profile.old)] = comp
  name.pos = setdiff(profile.pos, profile.old)
  name.neg = setdiff(profile.neg, profile.old)
  
  if(length(name.pos) != 0 & length(name.neg) != 0){
    var.sign[name.pos] = 'pos'
    var.sign[name.neg] = 'neg'
    mean.tmp = t(as.matrix(apply(data[,setdiff(profile, profile.old)],1,function(OTU)tapply(OTU, var.sign[setdiff(profile, profile.old)], mean))))
    colnames(mean.tmp) = paste0('comp', comp, colnames(mean.tmp)) 
    diablo.means = cbind(diablo.means, mean.tmp)
  }else if(length(name.pos) != 0){
    var.sign[name.pos] = 'pos'
    diablo.means = cbind(diablo.means,  apply(as.matrix(data[, name.pos]), 1, mean))
    colnames(diablo.means) = c(colnames(diablo.means)[-ncol(diablo.means)], paste0('comp', comp, 'pos'))
  }else{  # if length(name.pos) = 0
    var.sign[name.neg] = 'neg'
    diablo.means = cbind(diablo.means,  apply(as.matrix(data[, name.neg]), 1, mean))
    colnames(diablo.means) = c(colnames(diablo.means)[-ncol(diablo.means)], paste0('comp', comp, 'neg'))
  }
  profile.old = c(profile.old, profile)
}

# this shows the mean profile plots within the pos / neg profiles, per component
# careful the scale can be misleading!
col.remove = unique(which(is.na(diablo.means), arr.ind = TRUE)[, 'col'])
if(length(col.remove) >0) diablo.means = diablo.means[, -c(col.remove)]


# all profile plots, on the original variables, selected
# remove variables not selected
var.selected = var.selected[!is.na(var.selected)]
#length(var.selected)   
var.sign = var.sign[!is.na(var.sign)]
#length(var.sign)

return(list(var.selected = var.selected, var.sign = var.sign, diablo.means = diablo.means))

}
#############################################################################################################
# Authors:
#  Jasmin Straube, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2016
# last modified: 2016
#
# Copyright (C) 2016
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# 24 Oct, KA: debugging as original function was filtering on the samples, not the variables
# liver tocixity
# data = data.liver
# group = group.liver
# time = time.liver
# dim(meltd)     #[1] 199424      3

# TF data
# data = data.TF
# group = NULL
# time = time.TF
##==========================Filtering of low fold changes across time =========================##

#see remark below RE foldchange (minmax calc) not clear
#filter.foldchange <- function(data,time,group,threshold=log2(1.5),logFC=TRUE){
filter.foldchange <- function(data,time,group, threshold = 1, logFC=TRUE){
  data.temp <- data.frame(data)
  #if(!is.null(group)){  # KA added if condition when there is only one group
  data.temp$group <- paste(group,'/',time)
  #}else{
  #  data.temp$group <- time
  #}
  
  # melt data with respect to group * time
  meltd <- melt(data.temp,id=c("group"),na.rm=T)
  
  # calculate mean across group * time
  #if(!is.null(group)){  # KA added if condition when there is only one group
  cast <- dcast(data = meltd,group~variable,fun.aggregate = mean,na.rm=T)
  # }else{
  #   meltd <- melt(data.temp,id=c("time"),na.rm=T)
  #   cast <- dcast(data = meltd,group~time,fun.aggregate = mean,na.rm=T)
  # }
  
  # extract group ID? (flo to check: is there an easier coding than this?)
  if(!is.null(group)){  # KA added if condition when there is only one group
    id <- sapply(cast[,1],function(x)strsplit(as.character(x),split = '/')[[1]][1])
  }else{
    id <- sapply(cast[,1],function(x)strsplit(as.character(x),split = '/')[[1]][2])
  }
  
  if(!is.null(group)){ # KA added if condition when there is only one group
  if(logFC){  # if data have been log transformed
    # calculate the log FC between the min and the max mean values for each variable, and each time (?)
    # (I dont understand that function - KA)
  minmax <- apply(cast[,-1],2,function(x)tapply(x,id,function(y)max(y,na.rm = T)-min(y,na.rm = T)))
  }else{ # for the case of no group, but time points, not sure this is doing what it is supposed to do:
    # in fact what happens if there is a zero value here?
    minmax <- apply(cast[,-1],2,function(x)tapply(x,id,function(y)max(y,na.rm = T)/min(y,na.rm = T))) 
  }
  }else{ # only one group
    if(logFC){  # if data have been log transformed
      minmax <- apply(cast[,-1],2,function(x){max(x,na.rm = T)-min(x,na.rm = T)})
    }else{ # for the case of no group, but time points, not sure this is doing what it is supposed to do:
      # in fact what happens if there is a zero value here?
      minmax <- apply(cast[,-1],2,function(x){max(x,na.rm = T)/min(x,na.rm = T)})
    }
    
  }

  #need justification for filtering
  #could use a quantile
  #q <- quantile(as.vector(minmax),0.40)
  
  #remove the variables that fall below the defined fc threshold for ALL groups
  # !! KA rk: I think that should be threshold = 1.5 instead of log2(1.5)?
  if(!is.null(group)){  # KA added if condition when there is only one group
    index <- which(colSums(minmax<threshold)==length(unique(group)))
  }else{
    index <- which(minmax<threshold)
  }
  
  #l <- list(filtered.data = data[-index,],filter.index=index, foldChange=minmax, threshold=threshold,log=logFC)
  # the data are n x p, so we remove on the columns
  l <- list(filter.data = data[, -index], filter.index=index, foldChange=minmax, threshold=threshold,log=logFC)
  class(l) <- 'fc'
  return(l)
}

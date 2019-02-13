\name{filter.foldchange}
\encoding{latin1}
\alias{filter.foldchange}	

\title{Filtering of variables with low fold change}

\description{
The function removes the variables that fall below a pre-defined threshold.
}

\usage{
filter.foldchange(data,
    group,
    time,
	threshold=log2(1.5),
    logFC = TRUE)
}

\arguments{
\item{data}{ numeric \code{data.frame} or \code{matrix}, containing measurements of a multivariate time series. Sample information are in row and variables (e.g. genes) in columns.}
  \item{group}{ \code{vector} of length \code{nrow(data)} containing the grouping information of the samples (e.g. healthy, control).}
  \item{time}{ \code{vector} of length \code{nrow(data)} containing the time point information of the samples (e.g. 1,2,...). }
  \item{threshold}{ \code{numeric} value indicating the fold change threshold. Variables that fall below the threshold for all groups will be removed. Default value is log2(1.5)} % need to change that
  \item{logFC}{ \code{logical} value indicating if the data were log transformed. Default value is \code{TRUE}.}
  }

\details{
The purpose of this approach is to remove variables that show little variation across time points.

For each variable, time point and group the mean is extracted. When logFC = T the fold change of a variable is calculated as the difference 
of the maximum and minimum mean of the time points. If logFC = F the fold change is the ratio of maximum and minimum mean of the time points.

Variables are considered not informative if they do not pass the defined \code{threshold} across all groups. The index of these variables are reported in \code{filter.index}.
}

\value{
\code{filter.foldchange} returns an object of class \code{fc}, a list that 
contains the following components:
  \item{filter.data}{ \code{data.frame} the original data without the variables that did not pass the threshold.}
  \item{filter.index}{ \code{vector} containing the indices of the variables that where removed from the original data set.}
  \item{foldChange}{ a \code{vector} or \code{matrix} containing the maximum average (log) fold change across time for each variable (row) and group (column).}
  \item{threshold}{ a \code{numeric} value indicating the threshold applied to filter variables.}
  \item{logFC}{ a \code{logical} value indicating if the log fold change (\code{TRUE}) or fold change (\code{FALSE}) is returned.}
}

\author{Jasmin Straube.}

\seealso{\code{\link{plot.fc}}, \code{\link{lmmSpline}}}

\examples{
\dontrun{
#Example
data("liver.toxicity")
data <- liver.toxicity$gene
time <- liver.toxicity$treatment$Time.Group
group <- liver.toxicity$treatment$Dose.Group
filter.data <-  filter.foldchange(data,group=group,time=time,threshold=log2(1.1),logFC=T)$filtered.data
}}



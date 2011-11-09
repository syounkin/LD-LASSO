\name{block.bounds}
\alias{block.bounds}
\title{ Creates a vector of haplotype block boundary basepair positions. }
\description{
This function takes a map vector and a block coordinate vector and
creates a vector of block boundary positions for use when plotting block
boundaries, as in plot_ldlasso.
}
\usage{
block.bounds(map, block.cood)
}
\arguments{
  \item{map}{
    a vector of SNP map positions in basepairs (or kilo-basepairs)
  }
  \item{block.cood}{
  A vector of length p+1, where p is the number of SNPs.  block.cood is an indicator vector
  that indicates block boundaries at all p+1 SNP bounded intervals.  Use
  find.bounds to create this vector.
  }
}

\value{
A vector of block boundary positions in basepairs
}
\author{
Samuel G. Younkin
}
\examples{
data("ldlasso_example")
bpmap <- (block.obj@gtdata@map-block.obj@gtdata@map[1])/1e3
block.bounds.vec <- block.bounds( map = bpmap , block.cood = block.cood )
block.bounds.vec

\dontshow{

function( map, block.cood ){
  x.bd <- c()
  wd <- map[length(map)]-map[1]
  for( i in 1:(length(map)+1) ){
    if( i == 1 ){
      x.bd[i] <- map[1]
    }else if( any(i==2:length(map) ) ){
      x.bd[i] <- mean(map[(i-1):i])
    }else if( i==(length(map)+1)){
      x.bd[length(map)+1] <- map[length(map)]
    }
  }
  x.bd[-c(1,length(map)+1)][block.cood[-c(1,length(map)+1)]==1]
  }
}
}

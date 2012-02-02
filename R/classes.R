setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass( "ldlasso", representation(
                                    geno = "matrix",
                                    pheno = "integer",
                                    s1 = "numericOrNULL",
                                    s2 = "numeric",
                                    r2 = "numeric",
                                    beta = "numeric",
                                    OR = "numeric",
                                    logOR.norm = "numeric",
                                    var = "numeric",
                                    maf.case = "numeric",
                                    maf.con = "numeric",
                                    maf.tot = "numeric"
                                    )
         )

ldlasso =
# an ldlasso object with genotype matrix (subjects by SNPs) and phenotype vector
function( geno, pheno, s1, s2, r2 )
{
  geno = as( geno, "matrix" )
  pheno = as( pheno, "integer" )
  logOR.obj <- logOR( geno, pheno )
  if( nrow( geno) != length(pheno) )
    stop( "The number of rows in geno should equal the length of pheno" )
  if(is.null(colnames(geno)))
     colnames(geno) <- paste( "SNP", 1:ncol(geno), sep = "" )
  if(is.null(rownames(geno)))
     rownames(geno) <- paste( "sub", 1:nrow(geno), sep = "" )
  new( "ldlasso", geno = geno, pheno = pheno, s1 = s1, s2 = s2, r2 = r2, OR = logOR.obj$OR, logOR.norm = logOR.obj$y, var = logOR.obj$var_y, maf.case = logOR.obj$f2, maf.con = logOR.obj$f1, maf.tot = logOR.obj$f0 )
}



n.sub <- 1000; n.snp <- 50;
geno <- matrix( sample( 0:2, n.sub*n.snp, replace = TRUE ), nrow = n.sub, ncol = n.snp )
pheno <- sample( 0:1, n.sub, replace = TRUE )

ldlasso.obj <- ldlasso( geno = geno, pheno = pheno, s1 = NULL, s2 = 1, r2 = 0 )
ldlasso.obj <- findS1(ldlasso.obj)

ldlasso.obj

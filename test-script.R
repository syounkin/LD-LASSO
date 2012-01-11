library("ldlasso")
data("ldlasso-example")
ldlasso.obj <- ldlasso( geno = geno, pheno = pheno, s1 = 1, s2 = 1e-2, r2 = 0 )
ldlasso.obj@beta <- solve(ldlasso.obj)
ldlasso.obj

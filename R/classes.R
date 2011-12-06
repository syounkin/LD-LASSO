setClass("ldlasso", representation( geno = "matrix", pheno = "vector", s1 = "numeric", s2 = "numeric", r2 = "numeric", delta = "numeric" ), prototype = list( delta = 1e-6 )  )

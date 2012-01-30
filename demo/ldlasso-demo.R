## n.sub <- 1000; n.snp <- 50;
## geno <- matrix( sample( 0:2, n.sub*n.snp, replace = TRUE ), nrow = n.sub, ncol = n.snp )
## pheno <- sample( 0:1, n.sub, replace = TRUE )

## ldlasso.obj <- ldlasso( geno = geno, pheno = pheno, s1 = NULL, s2 = 1, r2 = 0 )
## ldlasso.obj <- findS1(ldlasso.obj)

## ldlasso.obj

library("CC.Sim")
data("sample")
geno.mat <- geno.mat[1:100,-17]

# Simulation parameters
causal.index <- 13
p0 <- 0.10 # probability of disease given homo major
beta <- 100 # "effect size"

hapslist <- GetHaps( geno.mat = geno.mat)
haps <- hapslist$haps
hap.probs <- hapslist$hap.probs
p <- sum((haps[,causal.index]==1)*hap.probs) # major allele frequency of causal SNP

problist <- GenotypeProb(p0 = p0, beta = beta, p = p)
pgd <- problist$pgd # probability of genotype, given disease
pgh <- problist$pgh # probability of genotype, given healthy

# Simulate Case-Control data set with 'size' of each
size <- 1e3
sim.data.list <- CC.Sim( pgd = pgd, pgh = pgh, haps = haps, hap.probs = hap.probs, causal.index = causal.index, size = size )

geno <- rbind( sim.data.list[[1]], sim.data.list[[2]] )
pheno <- c( rep(1, size), rep(0, size) )

ldlasso.obj <- ldlasso( geno = geno, pheno = pheno, s1 = NULL, s2 = 1e1, r2 = 0 )

findS1(ldlasso.obj)

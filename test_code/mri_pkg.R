# Mesoscale real data analysis (Auckland mri data)

#### load networks and metadata ####

# install package
library(devtools)
install_github('peterwmacd/nettest')
library(nettest)
# wd
setwd('~/packages/nettest/test_code/')
source('aux_functions.R')
# study
study <- 'taowu'
# parcellation + number of regions
parcels <- 'AAL116'
n <- 116

# load metadata
metadata <- read.csv(paste0('auckland_mri_data/Metadata_v4/',study,'_metadata.csv'))
# total sample size
m1 <- sum(metadata$Group=='Control')
m2 <- sum(metadata$Group=='PD')
M <- m1 + m2

# control networks
A1 <- array(NA,c(n,n,m1))
A1l <- list()
ids1 <- metadata$Subject[metadata$Group=='Control']
for(kk in 1:m1){
  matfile <- paste0('auckland_mri_data/',study,'/sub-',ids1[kk],
                    '/sub-',ids1[kk],'_AAL116_correlation_matrix.mat')
  temp <- R.matlab::readMat(matfile)
  A1[,,kk] <- A1l[[kk]] <- hollowone(temp$data)
}

# case networks
A2 <- array(NA,c(n,n,m2))
A2l <- list()
ids2 <- metadata$Subject[metadata$Group=='PD']
for(kk in 1:m2){
  matfile <- paste0('auckland_mri_data/',study,'/sub-',ids2[kk],
                    '/sub-',ids2[kk],'_AAL116_correlation_matrix.mat')
  temp <- R.matlab::readMat(matfile)
  A2[,,kk] <- A2l[[kk]] <- hollowone(temp$data)
}

# load parcellation metadata
# required package brainGraph
metadata_parcel <- brainGraph::aal116
# breakpoints for lobes
lbreak <- which(c(metadata_parcel$lobe,NA)!=c(NA,metadata_parcel$lobe))-.5

# indices for hypothesis blocks
hfl <- which(metadata_parcel$lobe=='Frontal')
hcb <- which(metadata_parcel$lobe=='Cerebellum')
hall <- c(hfl,hcb)

# rectangular specifications/masking sets for each test:

# test1: FL x FL
rect1 <- list(hfl,hfl)
mask1 <- list(list(hcb,hcb),list(hfl,hcb))

# test2: CB x CB
rect2 <- list(hcb,hcb)
mask2 <- list(list(hfl,hfl),list(hfl,hcb))

# test3: FL x CB
rect3 <- list(hfl,hcb)
mask3 <- list(list(hfl,hfl),list(hcb,hcb))

#### testing results ####

# for three blocks, store holm adjusted -log(p-values) (and sds for random projection)
# 1 basic approach
# 5 projection tests
# 5 random projection tests (100 reps)
dvec <- 2*1:5 # dimensions for mesoscale testing

pval_table <- matrix(NA,12,6)
rownames(pval_table) <- c('basic',
                          paste0('proj. d=',dvec),
                          paste0('rand. proj. d=', dvec),
                          'block')
colnames(pval_table) <- c('FL/FL','CB/CB','FL/CB',
                          'rej. FL/FL','rej. CB/CB','rej. FL/CB')

#### basic testing ####

# Basic test: ignores the network structure of the data, just looks at the
# edges in the block of interest as a vector of observations and tests for
# any differences in means... (with an F test)

# 1. diagonal block: FL <-> FL
bt1 <- basic_ftest(A1,A2,hfl,hfl,0.1)

# 2. diagonal block: CB <-> CB
bt2 <- basic_ftest(A1,A2,hcb,hcb,0.1)

# 3. off-diagonal block: FL <-> CB
bt3 <- basic_ftest(A1,A2,hfl,hcb,0.1)

# then get the adjusted fdr corrected p-values
pb_bon <- p.adjust(c(bt1$pval,bt2$pval,bt3$pval),method='none')
# same conclusion, only reject the CB <-> CB block hypothesis at level 0.1, same
# conclusion with a holm correction
# store results
pval_table[1,1:3] <- pb_bon
pval_table[1,4:6] <- as.integer(pb_bon < 0.05)

#### block testing ####

# Block testing: testing after a one dimensional projection by comparing the
# means of each block (note this is a special case of the mesoscale test with
# prespecified projections), this reduces dimension but 'naively', and may lose
# power relative to the basic test

Pfl <- matrix(0,n,1)
Pfl[hfl] <- 1/sqrt(length(hfl))
# sub projection (o/n basis) for cerebellum
Pcb <- matrix(0,n,1)
Pcb[hcb] <- 1/sqrt(length(hcb))

# 1. diagonal block: FL <-> FL
blt1 <- nettest:::Weight_meso(A1l,A2l,
                              hyp_indices=as.matrix(expand.grid(rect1)),
                              hyp_proj=list(Lproj=Pfl,Rproj=Pfl),
                              var_type='basic',
                              directed=FALSE)

# 2. diagonal block CB <-> CB
blt2 <- nettest:::Weight_meso(A1l,A2l,
                              hyp_indices=as.matrix(expand.grid(rect2)),
                              hyp_proj=list(Lproj=Pcb,Rproj=Pcb),
                              var_type='basic',
                              directed=FALSE)

# 3. off-diagonal block FL <-> CB
blt3 <- nettest:::Weight_meso(A1l,A2l,
                              hyp_indices=as.matrix(expand.grid(rect3)),
                              hyp_proj=list(Lproj=Pfl,Rproj=Pcb),
                              var_type='basic',
                              directed=FALSE)

# then get the adjusted fdr corrected p-values
pb_block <- p.adjust(c(blt1$pval,blt2$pval,blt3$pval),method='none')
# same conclusion, only reject the CB <-> CB block hypothesis at level 0.1, same
# conclusion with a holm correction
# store results
pval_table[12,1:3] <- pb_block
pval_table[12,4:6] <- as.integer(pb_block < 0.05)

#### projection testing ####

# Mesoscale testing, here we test the hypothesis for each block by learning a
# d-dimensional projection from the rest of the network, for d=2,4,6,8,10.

for(dd in dvec){
  resrow <- (dd/2)+1

  # 1. diagonal block FL <-> FL
  # pt1 <- proj_ftest(A1,A2,hfl,hfl,0.1,
  #                   Pfl,Pfl,tilde=FALSE)
  pt1 <- nettest::Mesoscale_test(A1l,A2l,sig=0.1,
                                 hyp_set=rect1,
                                 dimension=dd,
                                 directed=FALSE,self_loops=FALSE,
                                 proj_type='one_step',
                                 masked_set=mask1)

  # 2. diagonal block CB <-> CB
  pt2 <-  nettest::Mesoscale_test(A1l,A2l,sig=0.1,
                                  hyp_set=rect2,
                                  dimension=dd,
                                  directed=FALSE,self_loops=FALSE,
                                  proj_type='one_step',
                                  masked_set=mask2)

  # 3. off-diagonal block FL <-> CB
  pt3 <- nettest::Mesoscale_test(A1l,A2l,sig=0.1,
                                 hyp_set=rect3,
                                 dimension=dd,
                                 directed=FALSE,self_loops=FALSE,
                                 proj_type='one_step',
                                 masked_set=mask3)

  # p-value adjustment
  pp_bon <- p.adjust(c(pt1[2],pt2[2],pt3[2]),method='none')

  pval_table[resrow,1:3] <- pp_bon
  pval_table[resrow,4:6] <- as.integer(pp_bon < 0.05)
  print(paste0('done ',dd))
}

#### random projection testing ####

# Another comparison for mesoscale testing: do the same approach but specify a
# random projection of the data for d=2,4,6,8,10 for 50 replications.
# Shows that the mesoscale approach is indeed getting useful low-dimensional
# information out of the data.

randreps <- 50
pvals_rand <- array(NA,c(length(dvec),3,randreps))
set.seed(20)

for(dd in dvec){
  resrow <- (dd/2)

  for(bb in 1:randreps){
    # sub projection (o/n basis) for frontal lobe
    Prand <- svd(matrix(rnorm(n*dd),ncol=dd))$u

    # 1. diagonal block FL <-> FL
    pt1 <- nettest:::Weight_meso(A1l,A2l,
                                 hyp_indices=as.matrix(expand.grid(rect1)),
                                 hyp_proj=list(Lproj=Prand,Rproj=Prand),
                                 var_type='basic',
                                 directed=FALSE)

    # 2. diagonal block CB <-> CB
    pt2 <- nettest:::Weight_meso(A1l,A2l,
                                hyp_indices=as.matrix(expand.grid(rect2)),
                                hyp_proj=list(Lproj=Prand,Rproj=Prand),
                                var_type='basic',
                                directed=FALSE)

    # 3. off-diagonal block FL <-> CB
    pt3 <-  nettest:::Weight_meso(A1l,A2l,
                                  hyp_indices=as.matrix(expand.grid(rect3)),
                                  hyp_proj=list(Lproj=Prand,Rproj=Prand),
                                  var_type='basic',
                                  directed=FALSE)

    # p-value adjustment
    pp_bon <- p.adjust(c(pt1$pval,pt2$pval,pt3$pval),method='none')

    pvals_rand[resrow,,bb] <- pp_bon
  }

  for(ii in 1:3){
    # report median, and proportion of trials rejected at level 0.05
    pval_table[(resrow+6),ii] <- median(pvals_rand[resrow,ii,])
    pval_table[(resrow+6),(ii+3)] <- mean(pvals_rand[resrow,ii,] < 0.05)
  }

  print(paste0('done ',dd))
}
# takes a couple minutes to run

pval_table <- round(pval_table,7)
# report as a table in the manuscript

# save results
saveRDS(pval_table,file='mri_pkg_table.rds')

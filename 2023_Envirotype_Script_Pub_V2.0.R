# Envirotypic Genomic Prediction
# Noah D. Winans


rm(list = ls()) # clear global environment
# load("C:............/output_for_pub_Script.RData") # Load data file

#Units for plant height is given in inches, yield is given in lbs/ac
#Note that 18LY as presented in the manuscript is shown as 18RF in the data
#K_A, K_D, K_T, and K_W represent relationship matrices based on 
#additive, dominance, envirotypic T, and envirotypic W respectively.
#T_Matrix2 and W_matrix represent the envirotypic data used to create the
#envirotypic relationship matrices. The first ten columns are soil characteristics

# Load packages 
library(progress)
library(fastmatrix)
library(BGLR)
library(dplyr)


##### Build Kernels #####

# Envirotypic Kernels:
jp <- matrix(1, nrow = 100, ncol = 100)
Ze <- model.matrix(~ env-1, Phenotype_data2)
Za <- model.matrix(~ gid-1, Phenotype_data2)

colnames(Ze)
colnames(Ze) <- c("18COL", "18CS", "18GC", "18RF", "18VC", "19COL", "19CS", "19HAY", "19TA", "19VC") 

Ze <- Ze[, unique(Phenotype_data2$env)]

left_add <- Ze%*%K_W%*%t(Ze)
right_add <- Za%*%K_A%*%t(Za)
right_dom <- Za%*%K_D%*%t(Za)

KGE_AW <- hadamard(left_add, right_add)
KGE_DW <- hadamard(left_add, right_dom)
KE_W <- kronecker(K_W, jp)

left_add<-Ze%*%K_T%*%t(Ze)
right_add<-Za%*%K_A%*%t(Za)
right_dom<-Za%*%K_D%*%t(Za)

KGE_AT <- hadamard(left_add, right_add)
KGE_DT <- hadamard(left_add, right_dom)
KE_T <- kronecker(K_T, jp)

# Other Kernels
jq <- matrix(1, nrow = 10, ncol = 10)

left_add <- Ze%*%diag(1,10,10)%*%t(Ze)
right_add <- Za%*%K_A%*%t(Za)
right_dom <- Za%*%K_D%*%t(Za)

KG_GE_A <- hadamard(left_add,right_add)
KG_GE_D <- hadamard(left_add,right_dom)

KG_G_A <- kronecker(jq,K_A)
KG_G_D <- kronecker(jq,K_D)

#### Building Models for BGLR ####

Ze <- Ze[, unique(Phenotype_data2$env)]

M5_T <- list(list(K = KG_G_A, model = "RKHS"),
             list(K = KG_G_D, model = "RKHS"),
             list(K = KE_T, model = "RKHS"),
             list(K = KGE_AT, model = "RKHS"),
             list(K = KGE_DT, model = "RKHS"),
             list(X = Ze, model = 'BRR'))

M4_T <- list(list(K = KG_G_A, model = "RKHS"),
             list(K = KG_G_D, model = "RKHS"),
             list(K = KE_T, model = "RKHS"),
             list(X = Ze, model = 'BRR'))

M5_W <- list(list(K =KG_G_A, model = "RKHS"),
             list(K = KG_G_D, model = "RKHS"),
             list(K = KE_W, model = "RKHS"),
             list(K = KGE_AW, model = "RKHS"),
             list(K = KGE_DW, model = "RKHS"),
             list(X = Ze, model = 'BRR'))

M4_W <- list(list(K = KG_G_A, model = "RKHS"),
             list(K = KG_G_D, model = "RKHS"),
             list(K = KE_W, model = "RKHS"),
             list(X = Ze, model = 'BRR'))

M3 <- list(list(K = KG_G_A, model = "RKHS"),
           list(K = KG_G_D, model = "RKHS"),
           list(K = KG_GE_A, model = "RKHS"),
           list(K = KG_GE_D, model = "RKHS"),
           list(X = Ze, model = 'BRR'))

M2 <- list(list(K = KG_G_A, model = "RKHS"),
           list(K = KG_G_D, model = "RKHS"),
           list(X = Ze, model = 'BRR'))

M1 <- list(list(K = KG_G_A, model = "RKHS"),
           list(X = Ze, model = 'BRR'))


#### Cross Validation Schemes ####

# Replace Value in Phenotype_data2 with trait of interest from data2
Phenotype_data2$value <- data2$Yield

######  CV0: Predicting Unobserved Environments ######
environments <- c(unique(Phenotype_data2$env), unique(Phenotype_data2$env))
set.seed(123)
cycles <- 20
ETAX <- M2 #Specify Model Here
pb <- progress_bar$new(total = cycles)
CV0.Acc.Out <- matrix(nrow = cycles, ncol = 1)

for(r in 1:cycles){
  test_env <- environments[r]
  train_env <- setdiff(environments,test_env)
  Phenotype_CV0 <- Phenotype_data2
  Phenotype_CV0$Y2 <- NA
  Phenotype_CV0$Y2[Phenotype_CV0$env%in%train_env] <- Phenotype_CV0$value[Phenotype_CV0$env%in%train_env]
  y_t <- as.numeric(Phenotype_CV0$Y2)
  fit <- BGLR(y = y_t, ETA = ETAX, nIter = 10000, burnIn = 1000, verbose = F) # Fit Model
  real <- as.numeric(Phenotype_CV0$value[fit$whichNa])
  pred <- fit$yHat[fit$whichNa]
  CV0.Acc.Out[r, 1] <- cor(pred, real, use = "complete.obs")
  pb$tick()
}

  ###### CV1: Predicting New Hybrids ######

set.seed(123)
cycles<-25
hybrid<-unique(Phenotype_data2$gid)
ETAX<-M1 #Specify Model Here
Z<-30 #Specify how many hybrids in validation
pb <- progress_bar$new(total = cycles)
CV1.Acc.Out <- matrix(nrow = cycles, ncol=1)
for(r in 1:cycles){
  test_geno <- sample(hybrid,Z)
  train_geno <- setdiff(hybrid,test_geno)
  Phenotype_CV1 <- Phenotype_data2
  Phenotype_CV1$Y2 <- NA
  Phenotype_CV1$Y2[Phenotype_CV1$gid%in%train_geno] <- Phenotype_CV1$value[Phenotype_CV1$gid%in%train_geno]
  y_t <- as.numeric(Phenotype_CV1$Y2)
  fit <- BGLR(y = y_t, ETA = ETAX, nIter = 10000, burnIn = 1000, verbose = F)
  real <- as.numeric(Phenotype_CV1$value[fit$whichNa])
  pred <- fit$yHat[fit$whichNa]
  CV1.Acc.Out[r, 1] <- cor(pred, real, use = "complete.obs")
  pb$tick()
}

###### CV2: Predicting Sparse Hybrid Trials ######

set.seed(123)
env <- as.numeric(length(levels(as.factor(Phenotype_data2$env))))
cycles <- 25
X <- 500  #Number in test sample
Z <- M1 #Which model
pb <- progress_bar$new(total = cycles)
CV2.Acc.Out <- matrix(nrow = cycles, ncol = 1)

for(r in 1:cycles){
  repeat{
    test_sample <- sample(1:1000,X)
    train_sample <- setdiff(1:1000,test_sample)
    Phenotype_CV2 <- Phenotype_data2
    Phenotype_CV2$Y2 <- NA
    Phenotype_CV2$Y2[train_sample] <- Phenotype_CV2$value[train_sample]
    na_genos<-Phenotype_CV2$gid[test_sample]
    if(all(count(as.data.frame(na_genos), na_genos) != env)){
      break
    }
  }
  y_t <- as.numeric(Phenotype_CV2$Y2)
  fit <- BGLR(y = y_t, ETA = Z, nIter = 10000, burnIn = 1000, verbose = F)
  real <- as.numeric(Phenotype_CV2$value[fit$whichNa])
  pred <- fit$yHat[fit$whichNa]
  CV2.Acc.Out[r,1] <- cor(pred, real, use = "complete.obs")
  pb$tick()
} 


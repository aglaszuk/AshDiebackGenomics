# R scripts used to perform analyses

Short description:

1_PCA.r, performs PCA and visualize geographic distribution of mother trees
2_Structure_tess3.r, estimates spatial population structure using the TESS3 method
3_LinkageDisequilibrium.r, obtains average LD, recombination rate under drift equilibrium, and half decay distance from TASSEL output, and plots LD decay
4_TestPredictors.r, tests for associations among predictors used in the GEA models
5_GWAS_binomial.r, runs logistic regression using adegenet and glmnet when ADB damage is coded as a 0/1 trait
6_GWAS_multinomial_ADB.r, runs multinomial logistic regression using adegenet and glmnet when ADB damage is coded as a ordinal categorical variable with 1 to 6 categories
7_GWAS_multinomial_flushing.r, runs multinomial logistic regression using adegenet and glmnet for first and second flushing coded as a ordinal categorical variable with 1 to 6 categories
8_GEA.r, perfomrs latent factor mixed model and redundancy analysis to test for associations between genotype and climatic variables of interest

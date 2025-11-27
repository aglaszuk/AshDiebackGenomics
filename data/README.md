# Data description

- The "1_raw_data" directory contains a pdf file describing raw data statistics reported by ThermoFischer and the link to the dryad archive containing the data
- The "2_Axiom_qualityFilterDat" directory contains the genotype matrix after applying filtering trshholds reccomended by ThermoFisher
- The "3_ControlsRemoval_Maf_miss_mono_Filters" directory contains the genotype matrix (both VCF and TXT format) after applying MAF and missingness filters and removing monomorphic loci
- The "4_NewRefAnn_UniqueHits_NoHybridInd_GWAS+LD" directory contains the genotype matrix (both VCF and TXT format) after retaining SNPs that blast to the chromosome-level reference genome. This dataset was used for GWAS and GEA analyses and LD estimates
- The "5_PhenotypeClimDat" directory contains metadata including ADB damage assessments and climate data obtained from WorldClim

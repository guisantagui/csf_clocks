################################################################################
# Merge Bader, adni_emory, cruchaga emory, emory diversity and emory 3
# datasets.
################################################################################

if (!require("sva", quietly = TRUE)){
        BiocManager::install("sva", update = FALSE)
}
library(sva)
if (!require("impute", quietly = TRUE)){
        BiocManager::install("impute", update = FALSE)
}
library(impute)

if (!require("MsCoreUtils", quietly = TRUE)){
        BiocManager::install("MsCoreUtils", update = FALSE)
}
library(MsCoreUtils)

if (!require("DreamAI", quietly = TRUE)){
        devtools::install_github("WangLab-MSSM/DreamAI/Code")
}

if(!require("glmnet", quietly = TRUE)){
        install.packages("glmnet")
}

library(DreamAI)
library(glmnet)

library(stringr)

# Functions
################################################################################

# Plot mean abundance vs missing proportion.
plot_na_distr <- function(prot_df, axis = 1){
        if (axis != 1){
                prot_df <- t(prot_df)
        }
        protein_means <- colMeans(prot_df, na.rm = TRUE)
        na_cnts <- colSums(is.na(prot_df), na.rm = TRUE)
        missing_pct <- na_cnts / nrow(prot_df)
        na_df <- data.frame(
                protein = colnames(prot_df),
                protein_means = protein_means,
                missing_count = na_cnts,
                missing_percentage = missing_pct)
        plot(protein_means, missing_pct,
             xlab = "Mean Log2 Intensity (When Detected)", 
             ylab = "Fraction of Missing Values",
             main = "MAR vs MNAR Diagnostic Plot")
        return(na_df)
}

get_best_alpha <- function(X, y, n_folds = 5, family = "gaussian"){
        alphas <- seq(0, 1, by = 0.1)
        
        best_alpha <- NULL
        best_score <- Inf
        
        for (a in alphas) {
                
                cvfit_tmp <- cv.glmnet(
                        x = as.matrix(X),
                        y = y,
                        family = family,
                        alpha = a,
                        nfolds = n_folds
                )
                
                score <- min(cvfit_tmp$cvm)   # inner CV MSE
                
                if (score < best_score) {
                        best_score <- score
                        best_alpha <- a
                }
        }
        out <- list(best_score = best_score,
                    best_alpha = best_alpha)
        return(out)
}

# File stuff
################################################################################

prot_info_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/timsTOF_filt_protInfo.tsv"
DGE_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/scRNAseq/DGE"


bader_prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/bader_prot_parsed.csv"
emoradni_prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/emory_csf_parsed.csv"
cruchaga_prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/cruchaga_csf_parsed.csv"
emory_div_prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/emory_diversity_prot_parsed.csv"
emory_3_prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/emory_3_parsed.csv"

bader_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/bader_mdat_parsed.csv"
emoradni_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/emory_samp_info.csv"
cruchaga_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/cruchaga_samp_info.csv"
emory_div_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/emory_diversity_mdat_parsed.csv"
emory_3_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/emory_3_mdat_parsed.csv"

# Load proteomics data and metadata
################################################################################

prot_info <- read.table(prot_info_f)

bader_prot <- plotUtils::read_table_fast(bader_prot_f, row.names = 1)
emoradni_prot <- plotUtils::read_table_fast(emoradni_prot_f, row.names = 1)
cruchaga_prot <- plotUtils::read_table_fast(cruchaga_prot_f, row.names = 1)
emory_div_prot <- plotUtils::read_table_fast(emory_div_prot_f, row.names = 1)
emory_3_prot <- plotUtils::read_table_fast(emory_3_prot_f, row.names = 1)

bader_mdat <- read.csv(bader_mdat_f)
emoradni_mdat <- read.csv(emoradni_mdat_f)
cruchaga_mdat <- read.csv(cruchaga_mdat_f)
emory_div_mdat <- read.csv(emory_div_mdat_f)
emory_3_mdat <- read.csv(emory_3_mdat_f)

head(bader_mdat)

# Load single cell data
################################################################################
cell_files <- list.files(DGE_dir)
cell_files <- cell_files[grep(".tsv.gz", cell_files)]

deg_list <- list()
for (f in cell_files){
        cell_type <- strsplit(f, split = '_')[[1]][2]
        degs <- read.table(file.path(DGE_dir, f), header = TRUE)
        degs$cell_type <- cell_type
        deg_list[[cell_type]] <- degs
}


log2FC_thrshld <- 2

deg_list <- do.call(rbind.data.frame, deg_list)
rownames(deg_list) <- 1:nrow(deg_list)
hist(deg_list$log2FoldChange)
deg_list <- deg_list[deg_list$log2FoldChange >= log2FC_thrshld, ]

deg_list <- deg_list[!deg_list$gene %in% deg_list$gene[duplicated(deg_list$gene)], ]

# Harmonize metadata
################################################################################

bader_mdat_parsed <- bader_mdat[, c("Run_timsTOFMeasurement", "ID", "Sex", "Age", "Diagnosis_group", "Sample_preparation_batch")]
colnames(bader_mdat_parsed) <- c("sampleID", "individualID", "sex", "age", "diagnosis", "internal_batch")
bader_mdat_parsed$race <- NA
bader_mdat_parsed$apoe <- NA
bader_mdat_parsed$substudy <- "Bader"
bader_mdat_parsed$diagnosis <- gsub(
        "Neurological Control", "Control", bader_mdat_parsed$diagnosis
)
bader_mdat_parsed$sampleID <- make.names(bader_mdat_parsed$sampleID)

emoradni_mdat_parsed <- emoradni_mdat[, c("ProteomicID", "RID", "sex", "age", "DIAGNOSIS", "DXPARK", "DXAD", "DXDEP", "DXMDUE", "orig_prot", "race", "ethnicity", "apoe")]

table(emoradni_mdat_parsed$DIAGNOSIS, emoradni_mdat_parsed$DXAD)
table(emoradni_mdat_parsed$DIAGNOSIS, emoradni_mdat_parsed$DXPARK)

emoradni_mdat_parsed$AD <- grepl("Alzheimer", emoradni_mdat_parsed$DXMDUE)

emoradni_mdat_parsed[emoradni_mdat_parsed$DXAD == "Yes" & !is.na(emoradni_mdat_parsed$DXAD), ]
emoradni_mdat_parsed$AD[emoradni_mdat_parsed$DXAD == "Yes" & !is.na(emoradni_mdat_parsed$DXAD)] <- TRUE
emoradni_mdat_parsed$AD[is.na(emoradni_mdat_parsed$DXMDUE) & is.na(emoradni_mdat_parsed$DXAD)] <- NA
emoradni_mdat_parsed$PD[emoradni_mdat_parsed$DXPARK == "Yes" & !is.na(emoradni_mdat_parsed$DXPARK)] <- TRUE

emoradni_mdat_parsed$diagnosis_parsed <- NA
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "CN"] <- "Control"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "MCI" & (emoradni_mdat_parsed$AD & (!emoradni_mdat_parsed$PD | is.na(emoradni_mdat_parsed$PD)))] <- "MCI_AD"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "MCI" & (emoradni_mdat_parsed$PD & (!emoradni_mdat_parsed$AD | is.na(emoradni_mdat_parsed$AD)))] <- "MCI_PD"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "MCI" & ((!emoradni_mdat_parsed$PD | is.na(emoradni_mdat_parsed$PD)) & (!emoradni_mdat_parsed$AD | is.na(emoradni_mdat_parsed$AD)))] <- "MCI_other"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "MCI" & ((emoradni_mdat_parsed$PD & !is.na(emoradni_mdat_parsed$PD)) & (emoradni_mdat_parsed$AD & !is.na(emoradni_mdat_parsed$AD)))] <- "MCI_mixed"

emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "Dementia" & (emoradni_mdat_parsed$AD & (!emoradni_mdat_parsed$PD | is.na(emoradni_mdat_parsed$PD)))] <- "Dementia_AD"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "Dementia" & (emoradni_mdat_parsed$PD & (!emoradni_mdat_parsed$AD | is.na(emoradni_mdat_parsed$AD)))] <- "Dementia_PD"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "Dementia" & ((!emoradni_mdat_parsed$PD | is.na(emoradni_mdat_parsed$PD)) & (!emoradni_mdat_parsed$AD | is.na(emoradni_mdat_parsed$AD)))] <- "Dementia_other"
emoradni_mdat_parsed$diagnosis_parsed[emoradni_mdat_parsed$DIAGNOSIS == "Dementia" & ((emoradni_mdat_parsed$PD & !is.na(emoradni_mdat_parsed$PD)) & (emoradni_mdat_parsed$AD & !is.na(emoradni_mdat_parsed$AD)))] <- "Dementia_mixed"

table(emoradni_mdat_parsed$diagnosis_parsed)



emoradni_mdat_parsed$race[
        emoradni_mdat_parsed$ethnicity == "Hispanic or Latino"
] <- "Hispanic"
emoradni_mdat_parsed$race <- gsub(
        "Black or African American", "Black", emoradni_mdat_parsed$race
)

emoradni_mdat_parsed <- emoradni_mdat_parsed[, c("ProteomicID", "RID", "sex", "age", "diagnosis_parsed", "orig_prot", "race", "apoe")]

colnames(emoradni_mdat_parsed) <- c("sampleID", "individualID", "sex", "age", "diagnosis", "internal_batch", "race", "apoe")
emoradni_mdat_parsed$substudy <- "ADNI_Emory"
emoradni_mdat_parsed$sex <- ifelse(emoradni_mdat_parsed$sex == "Female", "F", "M")


cruchaga_mdat_parsed <- cruchaga_mdat[, c("GUSPECID", "RID", "sex", "age", "DIAGNOSIS", "DXPARK", "DXAD", "DXDEP", "DXMDUE", "orig_prot", "race", "ethnicity", "apoe")]

table(cruchaga_mdat_parsed$DIAGNOSIS, cruchaga_mdat_parsed$DXAD)
table(cruchaga_mdat_parsed$DIAGNOSIS, cruchaga_mdat_parsed$DXPARK)

cruchaga_mdat_parsed$AD <- grepl("Alzheimer", cruchaga_mdat_parsed$DXMDUE)

cruchaga_mdat_parsed[cruchaga_mdat_parsed$DXAD == "Yes" & !is.na(cruchaga_mdat_parsed$DXAD), ]
cruchaga_mdat_parsed$AD[cruchaga_mdat_parsed$DXAD == "Yes" & !is.na(cruchaga_mdat_parsed$DXAD)] <- TRUE
cruchaga_mdat_parsed$AD[is.na(cruchaga_mdat_parsed$DXMDUE) & is.na(cruchaga_mdat_parsed$DXAD)] <- NA
cruchaga_mdat_parsed$PD[cruchaga_mdat_parsed$DXPARK == "Yes" & !is.na(cruchaga_mdat_parsed$DXPARK)] <- TRUE

cruchaga_mdat_parsed$diagnosis_parsed <- NA
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "CN"] <- "Control"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "MCI" & (cruchaga_mdat_parsed$AD & (!cruchaga_mdat_parsed$PD | is.na(cruchaga_mdat_parsed$PD)))] <- "MCI_AD"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "MCI" & (cruchaga_mdat_parsed$PD & (!cruchaga_mdat_parsed$AD | is.na(cruchaga_mdat_parsed$AD)))] <- "MCI_PD"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "MCI" & ((!cruchaga_mdat_parsed$PD | is.na(cruchaga_mdat_parsed$PD)) & (!cruchaga_mdat_parsed$AD | is.na(cruchaga_mdat_parsed$AD)))] <- "MCI_other"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "MCI" & ((cruchaga_mdat_parsed$PD & !is.na(cruchaga_mdat_parsed$PD)) & (cruchaga_mdat_parsed$AD & !is.na(cruchaga_mdat_parsed$AD)))] <- "MCI_mixed"

cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "Dementia" & (cruchaga_mdat_parsed$AD & (!cruchaga_mdat_parsed$PD | is.na(cruchaga_mdat_parsed$PD)))] <- "Dementia_AD"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "Dementia" & (cruchaga_mdat_parsed$PD & (!cruchaga_mdat_parsed$AD | is.na(cruchaga_mdat_parsed$AD)))] <- "Dementia_PD"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "Dementia" & ((!cruchaga_mdat_parsed$PD | is.na(cruchaga_mdat_parsed$PD)) & (!cruchaga_mdat_parsed$AD | is.na(cruchaga_mdat_parsed$AD)))] <- "Dementia_other"
cruchaga_mdat_parsed$diagnosis_parsed[cruchaga_mdat_parsed$DIAGNOSIS == "Dementia" & ((cruchaga_mdat_parsed$PD & !is.na(cruchaga_mdat_parsed$PD)) & (cruchaga_mdat_parsed$AD & !is.na(cruchaga_mdat_parsed$AD)))] <- "Dementia_mixed"

table(cruchaga_mdat_parsed$diagnosis_parsed)

cruchaga_mdat_parsed$race[
        cruchaga_mdat_parsed$ethnicity == "Hispanic or Latino"
] <- "Hispanic"
cruchaga_mdat_parsed$race <- gsub(
        "Black or African American", "Black", cruchaga_mdat_parsed$race
)

cruchaga_mdat_parsed <- cruchaga_mdat_parsed[, c("GUSPECID", "RID", "sex", "age", "diagnosis_parsed", "orig_prot", "race", "apoe")]

colnames(cruchaga_mdat_parsed) <- c("sampleID", "individualID", "sex", "age", "diagnosis", "internal_batch", "race", "apoe")
cruchaga_mdat_parsed$substudy <- "ADNI_Cruchaga"
cruchaga_mdat_parsed$sex <- ifelse(cruchaga_mdat_parsed$sex == "Female", "F", "M")

head(emory_div_mdat)
emory_div_mdat_parsed <- emory_div_mdat[, c("TMT.ID", "GUID", "Sex", "Age", "Diagnosis", "Study", "Race", "APOE")]
colnames(emory_div_mdat_parsed) <- c("sampleID", "individualID", "sex", "age", "diagnosis", "internal_batch", "race", "apoe")
emory_div_mdat_parsed$sex <- ifelse(emory_div_mdat_parsed$sex == "Female", "F", "M")
emory_div_mdat_parsed$diagnosis <- str_trim(emory_div_mdat_parsed$diagnosis)
emory_div_mdat_parsed$diagnosis[emory_div_mdat_parsed$diagnosis == ""] <- NA
table(emory_div_mdat_parsed$race)
emory_div_mdat_parsed$race[grepl("Black", emory_div_mdat_parsed$race)] <- "Black"
emory_div_mdat_parsed$race[grepl("White", emory_div_mdat_parsed$race)] <- "White"
emory_div_mdat_parsed$race[emory_div_mdat_parsed$race == ""] <- NA
emory_div_mdat_parsed$substudy <- "Emory_Diversity"

head(emory_3_mdat)
emory_3_mdat_parsed <- emory_3_mdat[, c("GUID", "GUID", "Sex", "Age", "Diagnosis", "CSFneat.Batch", "Race", "APOE")]
colnames(emory_3_mdat_parsed) <- c("sampleID", "individualID", "sex", "age", "diagnosis", "internal_batch", "race", "apoe")
emory_3_mdat_parsed$sampleID <- make.names(emory_3_mdat_parsed$sampleID)
emory_3_mdat_parsed$diagnosis <- gsub("Normal", "Control", emory_3_mdat_parsed$diagnosis)
emory_3_mdat_parsed$race <- ifelse(grepl("Black", emory_3_mdat_parsed$race), "Black", "White")
emory_3_mdat_parsed$sex <- ifelse(emory_3_mdat_parsed$sex == 0, "F", "M")

emory_3_mdat_parsed$substudy <- "Emory_3"

intersection_emory <- intersect(emory_3_mdat_parsed$individualID, emory_div_mdat_parsed$individualID)

emory_3_mdat_parsed[match(intersection_emory, emory_3_mdat_parsed$individualID), ]
emory_div_mdat_parsed[match(intersection_emory, emory_div_mdat_parsed$individualID), ]

# There's some overlap of individuals between emory 3 and emory diversity.
# Remove the overlapping ones from Emory 3

emory_3_mdat_parsed <- emory_3_mdat_parsed[!emory_3_mdat_parsed$individualID %in% emory_div_mdat_parsed$individualID, ]
emory_3_prot <- emory_3_prot[, match(emory_3_mdat_parsed$sampleID, colnames(emory_3_prot))]

plot(density(as.matrix(bader_prot), na.rm = TRUE))


# Process the datasets individually and merge
################################################################################

bader_prot_proc <- t(bader_prot[, bader_mdat_parsed$diagnosis == "Control"])

keep_thrshld <- .3
bader_prot_proc <- bader_prot_proc[,
        apply(
                bader_prot_proc,
                2,
                function(x) sum(is.na(x))/length(x)
        ) < keep_thrshld
]


bader_mdat_parsed_ctrls <- bader_mdat_parsed[bader_mdat_parsed$diagnosis == "Control", ]
any(is.na(bader_mdat_parsed_ctrls$age))

#sum(apply(bader_prot_proc[, common_prots], 2, function(x) sum(is.na(x))/length(x)) > .2)/length(common_prots)


plot(density(bader_prot_proc, na.rm = TRUE))

bader_prot_proc <- log2(bader_prot_proc)
plot(density(bader_prot_proc, na.rm = TRUE))

na_mask <- is.na(bader_prot_proc)

bader_prot_proc_naDF <- plot_na_distr(bader_prot_proc)

plot(
        colMeans(apply(bader_prot_proc, 2, function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)), na.rm = TRUE),
        bader_prot_proc_naDF$missing_percentage
)

set.seed(123)

bader_prot_proc <- t(impute.knn(t(bader_prot_proc))$data)

#bader_prot_proc <- t(
#        MsCoreUtils::impute_matrix(
#                t(bader_prot_proc),
#                method = "mixed",
#                randna = bader_prot_proc_naDF$protein_means >= mean(bader_prot_proc, na.rm = TRUE),
#                mar = "knn",
#                mnar = "MLE"
#        )
#)

plot(density(bader_prot_proc[na_mask]), main="Distribution of Imputed Values Only")

#boxplot(bader_prot_proc[, 642][na_mask[, 642]],
#        bader_prot_proc[, 642][!na_mask[, 642]])

#apply(bader_prot_proc, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE) - 1.5 * sd(x, na.rm = TRUE), x))
#min(bader_prot_proc[, 1], na.rm = TRUE)
#mean(bader_prot_proc[, 1], na.rm = TRUE) - 1.5 * sd(bader_prot_proc[, 1], na.rm = TRUE)

bader_prot_proc_naDF$mean_after_imp <- colMeans(bader_prot_proc)

plot(bader_prot_proc_naDF$mean_after_imp, bader_prot_proc_naDF$missing_percentage)


bader_prot_proc <- apply(
        bader_prot_proc,
        2,
        function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
)
plot(density(bader_prot_proc, na.rm = TRUE))
plot(colMeans(bader_prot_proc), bader_prot_proc_naDF$missing_percentage)



emoradni_prot[emoradni_prot == -4] <- NA
emoradni_prot_proc <- as.matrix(emoradni_prot[emoradni_mdat_parsed$diagnosis == "Control", ])
emoradni_mdat_parsed_ctrls <- emoradni_mdat_parsed[emoradni_mdat_parsed$diagnosis == "Control", ]
sum(apply(emoradni_prot_proc, 2, function(x) sum(is.na(x))/length(x)) < keep_thrshld)
emoradni_prot_proc <- emoradni_prot_proc[, apply(emoradni_prot_proc, 2, function(x) sum(is.na(x))/length(x)) < keep_thrshld]

any(is.na(emoradni_mdat_parsed_ctrls$age))

plot(density(emoradni_prot_proc, na.rm = TRUE))

hist(emoradni_prot_proc)
unique(emoradni_prot_proc[emoradni_prot_proc < 0])
min(emoradni_prot_proc, na.rm = TRUE)

emoradni_prot_proc <- log2(emoradni_prot_proc)
plot(density(emoradni_prot_proc, na.rm = TRUE))

emoradni_prot_proc_naDF <- plot_na_distr(emoradni_prot_proc)
emoradni_prot_proc <- t(impute.knn(t(emoradni_prot_proc))$data)
emoradni_prot_proc_naDF$mean_after_imp <- colMeans(emoradni_prot_proc)

plot(emoradni_prot_proc_naDF$mean_after_imp, emoradni_prot_proc_naDF$missing_percentage)


emoradni_prot_proc <- apply(
        emoradni_prot_proc,
        2,
        function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
)
plot(density(emoradni_prot_proc, na.rm = TRUE))
plot(colMeans(emoradni_prot_proc), emoradni_prot_proc_naDF$missing_percentage)


cruchaga_prot_proc <- as.matrix(cruchaga_prot[cruchaga_mdat_parsed$diagnosis == "Control", ])
cruchaga_mdat_parsed_ctrls <- cruchaga_mdat_parsed[cruchaga_mdat_parsed$diagnosis == "Control", ]
any(is.na(cruchaga_mdat_parsed_ctrls$age))
cruchaga_prot_proc <- cruchaga_prot_proc[, apply(cruchaga_prot_proc, 2, function(x) sum(is.na(x))/length(x)) < keep_thrshld]


plot(density(cruchaga_prot_proc, na.rm = TRUE))

cruchaga_prot_proc <- log2(cruchaga_prot_proc)
plot(density(cruchaga_prot_proc, na.rm = TRUE))
cruchaga_prot_proc <- t(impute.knn(t(cruchaga_prot_proc))$data)
cruchaga_prot_proc <- apply(
        cruchaga_prot_proc, 2, function(x) (x - mean(x))/sd(x)
)
plot(density(cruchaga_prot_proc))

emory_div_prot_proc <- t(emory_div_prot[, emory_div_mdat_parsed$diagnosis == "Control" & !is.na(emory_div_mdat_parsed$diagnosis)])
emory_div_mdat_parsed_ctrls <- emory_div_mdat_parsed[emory_div_mdat_parsed$diagnosis == "Control" & !is.na(emory_div_mdat_parsed$diagnosis), ]
any(is.na(emory_div_mdat_parsed_ctrls$age))
emory_div_prot_proc <- emory_div_prot_proc[, apply(emory_div_prot_proc, 2, function(x) sum(is.na(x))/length(x)) < keep_thrshld]
plot(density(emory_div_prot_proc, na.rm = TRUE))
emory_div_prot_proc <- t(impute.knn(t(emory_div_prot_proc))$data)
emory_div_prot_proc <- apply(emory_div_prot_proc, 2, function(x) (x - mean(x))/sd(x))
plot(density(emory_div_prot_proc, na.rm = TRUE))


density(bader_mdat_parsed_ctrls$age)
density(emoradni_mdat_parsed_ctrls$age)
density(cruchaga_mdat_parsed_ctrls$age)
density(emory_div_mdat_parsed_ctrls$age)

common_prots <- colnames(bader_prot_proc) %>%
        intersect(colnames(emoradni_prot_proc)) %>%
        intersect(colnames(cruchaga_prot_proc)) %>%
        intersect(colnames(emory_div_prot_proc))

prot_info$protein
prot_info$protein_unique <- gsub("\\-.*", "", prot_info$protein)
prot_dups <- prot_info$protein_unique[duplicated(prot_info$protein_unique)]
prot_info[prot_info$protein_unique %in% prot_dups, ]
prot_info_simp <- dplyr::distinct(prot_info, ensembl_gene_id, protein_unique)
colnames(prot_info_simp) <- c("ensembl_gene_id", "protein")

prot_info_simp$cell_type <- deg_list$cell_type[match(prot_info_simp$ensembl_gene_id, deg_list$gene)]

table(prot_info_simp$cell_type[prot_info_simp$protein %in% common_prots])

#kruskal.test(
#        c(bader_mdat_parsed_ctrls$age,
#          emoradni_mdat_parsed_ctrls$age,
#          cruchaga_mdat_parsed_ctrls$age,
#          emory_div_mdat_parsed_ctrls$age),
#        g = c(bader_mdat_parsed_ctrls$substudy,
#              emoradni_mdat_parsed_ctrls$substudy,
#              cruchaga_mdat_parsed_ctrls$substudy,
#              emory_div_mdat_parsed_ctrls$substudy)
#)

mdat_all <- rbind.data.frame(
        bader_mdat_parsed_ctrls,
        emoradni_mdat_parsed_ctrls,
        cruchaga_mdat_parsed_ctrls,
        emory_div_mdat_parsed_ctrls
)

prot_all <- rbind.data.frame(
        bader_prot_proc[, common_prots],
        emoradni_prot_proc[, common_prots],
        cruchaga_prot_proc[, common_prots],
        emory_div_prot_proc[, common_prots]
)

prot_all_pca <- prcomp(prot_all, scale. = TRUE, center = TRUE)
mdat_all$sample <- mdat_all$sampleID
plotUtils::doPCAMultiPlot(prot_all_pca, samp_info = mdat_all, col = "substudy", nComps = 5)


bader_pca <- prcomp(prot_all[mdat_all$substudy == "Bader", ], scale. = TRUE, center = TRUE)
plotUtils::doPCAMultiPlot(
        bader_pca,
        samp_info = mdat_all[mdat_all$substudy == "Bader", ],
        col = "internal_batch",
        nComps = 5
)


kruskal.test(
        mdat_all$age,
        g = mdat_all$substudy
)


cell_types <- unique(deg_list$cell_type)

prot_info_simp <- prot_info_simp[prot_info_simp$protein %in% common_prots, ]

#cv_res_list <- list()
#final_mod_res <- list()
glm_list <- list()
for(ct in cell_types){
        print(sprintf("Fitting %s model...", ct))
        results_glm <- c()
        r2_glm_cv <- c()
        results_xgb <- c()
        r2_xgb_cv <- c()
        
        X_cell_train <- prot_all[
                mdat_all$substudy != "ADNI_Cruchaga",
                prot_info_simp$protein[prot_info_simp$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]
        ]
        X_cell_test <- prot_all[
                mdat_all$substudy == "ADNI_Cruchaga",
                prot_info_simp$protein[prot_info_simp$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]
        ]
        
        #X_cell_train <- prot_all[
        #        mdat_all$substudy == "Bader",
        #        prot_info_simp$protein[prot_info_simp$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]
        #]
        #X_cell_test <- prot_all[
        #        mdat_all$substudy == "ADNI_Emory",
        #        prot_info_simp$protein[prot_info_simp$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]
        #]
        
        y_train <- mdat_all$age[mdat_all$substudy != "ADNI_Cruchaga"]
        y_test <- mdat_all$age[mdat_all$substudy == "ADNI_Cruchaga"]
        
        #y_train <- mdat_all$age[mdat_all$substudy == "Bader"]
        #y_test <- mdat_all$age[mdat_all$substudy == "ADNI_Emory"]
        best_alpha <- get_best_alpha(X_cell_train, y_train)
        
        cvfit <- cv.glmnet(
                x = as.matrix(X_cell_train),
                y = y_train,
                family = "gaussian",
                alpha = best_alpha$best_alpha,
                nfolds = 5
        )
        pred_glm <- predict(cvfit, newx = as.matrix(X_cell_test), s = "lambda.1se")
        correlation <- cor(pred_glm[, 1], y_test)
        to_bind <- data.frame(cell_type = ct,
                              correlation = correlation)
        glm_list[[ct]] <- to_bind
}

do.call(rbind.data.frame, glm_list)


X_train <- prot_all[
        mdat_all$substudy != "ADNI_Cruchaga",
]
X_test <- prot_all[
        mdat_all$substudy == "ADNI_Cruchaga",
]

X_bader <- as.data.frame(prot_all[
        mdat_all$substudy == "Bader",
])
X_emory <- as.data.frame(prot_all[
        mdat_all$substudy == "ADNI_Emory",
])

X_cruchaga <- as.data.frame(prot_all[
        mdat_all$substudy == "ADNI_Cruchaga",
])

X_bader$age <- mdat_all$age[mdat_all$substudy == "Bader"]
X_emory$age <- mdat_all$age[mdat_all$substudy == "ADNI_Emory"]
X_cruchaga$age <- mdat_all$age[mdat_all$substudy == "ADNI_Cruchaga"]

lm_bader <- lm(age ~ ., X_bader)
lm_emory <- lm(age ~ ., X_emory)
lm_cruchaga <- lm(age ~ ., X_cruchaga)

cor(coef(lm_bader), coef(lm_emory), use = "complete.obs")
cor(coef(lm_bader), coef(lm_cruchaga), use = "complete.obs")

y_train <- mdat_all$age[mdat_all$substudy != "ADNI_Cruchaga"]
y_test <- mdat_all$age[mdat_all$substudy == "ADNI_Cruchaga"]
best_alpha <- get_best_alpha(X_train, y_train)

cvfit_all <- cv.glmnet(
        x = as.matrix(X_train),
        y = y_train,
        family = "gaussian",
        alpha = best_alpha$best_alpha,
        nfolds = 5
)
pred_glm_all <- predict(cvfit_all, newx = as.matrix(X_test), s = "lambda.1se")
correlation <- cor(pred_glm_all[, 1], y_test)


best_alpha <- get_best_alpha(X_bader[, colnames(X_bader) != "age"], X_bader$age)

cvfit_all <- cv.glmnet(
        x = as.matrix(X_bader[, colnames(X_bader) != "age"]),
        y = X_bader$age,
        family = "gaussian",
        alpha = best_alpha$best_alpha,
        nfolds = 5
)
bader_coefs <- coef(cvfit_all)
bader_coefs <- bader_coefs[, 1][bader_coefs[, 1] != 0]


lm_bader <- lm(age ~ ., X_bader[, c(names(bader_coefs)[names(bader_coefs) != "(Intercept)"], "age")])
lm_emory <- lm(age ~ ., X_emory[, c(names(bader_coefs)[names(bader_coefs) != "(Intercept)"], "age")])
lm_cruchaga <- lm(age ~ ., X_cruchaga[, c(names(bader_coefs)[names(bader_coefs) != "(Intercept)"], "age")])

coefs_bader <- coef(lm_bader)
coefs_bader <- coefs_bader[names(coefs_bader) != "(Intercept)"]

coefs_emory <- coef(lm_emory)
coefs_emory <- coefs_emory[names(coefs_emory) != "(Intercept)"]

coefs_cruchaga <- coef(lm_cruchaga)
coefs_cruchaga <- coefs_emory[names(coefs_cruchaga) != "(Intercept)"]

cor(coefs_bader, coefs_emory, use = "complete.obs")
cor(coefs_bader, coefs_cruchaga, use = "complete.obs")

plot(coefs_bader, coefs_emory)
plot(coefs_bader, coefs_cruchaga)
cor(coef(lm_bader), coef(lm_cruchaga), use = "complete.obs")

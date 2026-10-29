################################################################################
# CSF clocks: do a test for fitting the models                                 #
################################################################################
library(dplyr)

if (!require("lmerTest", quietly = TRUE)){
        install.packages("lmerTest")
}


library(caret)
if (!require("rsample", quietly = TRUE)){
        install.packages("rsample")
}
if (!require("xgboost", quietly = TRUE)){
        install.packages("xgboost")
}

data_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader"

DGE_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/scRNAseq/DGE"

prot_f <- "timsTOF_filt.tsv"
mdat_f <- "timsTOF_filt_mdat.tsv"
prot_info_f <- "timsTOF_filt_protInfo.tsv"
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/bader_models/"

plotUtils::create_dir_if_not(out_dir)

mdat_all <- as.data.frame(
        readxl::read_xlsx(
                "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/1-s2.0-S0092867426001042-mmc2.xlsx",
                sheet = 4
        )
)


# Functions
################################################################################
read_table_fast <- function(f, row.names = NULL, sep = "auto", sep2 = "auto"){
        df <- data.frame(fread(f, sep = sep, sep2 = sep2))
        if (!is.null(row.names)){
                rownames(df) <- df[, row.names]
                df <- df[, -row.names]
        }
        return(df)
}

# Load data
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

prot <- read.table(file.path(data_dir, prot_f))
mdat <- read.table(file.path(data_dir, mdat_f))
prot_info <- read.table(file.path(data_dir, prot_info_f))


table(deg_list$cell_type)

table(deg_list$cell_type[deg_list$gene %in% prot_info$ensembl_gene_id])

table(deg_list$cell_type[deg_list$gene %in% prot_info$ensembl_gene_id])

prot_info$cell_type <- deg_list$cell_type[match(prot_info$ensembl_gene_id, deg_list$gene)]
table(prot_info$cell_type)

# Do some plotting to assess bias of batches

ggplot(mdat, mapping = aes(x = Sample_preparation_batch, y = Age)) + geom_boxplot() + geom_jitter()

kruskal.test(mdat$Age, mdat$Sample_preparation_batch)

table(mdat$Sample_preparation_batch)

#fit_knn_mod 


library(rsample)

y <- mdat$Age
batch <- mdat$Sample_preparation_batch

X <- log2(t(as.matrix(prot)))


stopifnot(nrow(X) == length(y))


set.seed(1)
folds <- group_vfold_cv(
        mdat,
        group = Sample_preparation_batch,
        v = 7
)


batch_normalize <- function(train_x, train_batch, test_x, test_batch) {
        
        #train_x <- X_train
        #train_batch <- batch_train
        #test_x <- X_test
        #test_batch <- batch_test
        
        global_mean <- colMeans(train_x, na.rm = TRUE)
        global_sd   <- apply(train_x, 2, sd, na.rm = TRUE)
        
        split_idx <- split(seq_len(nrow(train_x)), train_batch)
        
        batch_means <- lapply(split_idx, function(idx) {
                colMeans(train_x[idx, , drop = FALSE], na.rm = TRUE)
        })
        
        batch_sds <- lapply(split_idx, function(idx) {
                apply(train_x[idx, , drop = FALSE], 2, sd, na.rm = TRUE)
        })
        
        norm_transform <- function(x, b) {
                mu_b <- batch_means[[b]]
                sd_b <- batch_sds[[b]]
                
                # handle NA / zero sd
                sd_b[sd_b == 0 | is.na(sd_b)] <- 1
                
                x2 <- sweep(x, 2, mu_b, "-")
                x2 <- sweep(x2, 2, sd_b, "/")
                x2 <- sweep(x2, 2, global_sd, "*")
                x2 <- sweep(x2, 2, global_mean, "+")
                x2
        }
        
        train_out <- do.call(rbind, lapply(unique(train_batch), function(b) {
                norm_transform(train_x[train_batch == b, , drop=FALSE], b)
        }))
        
        test_out <- do.call(rbind, lapply(unique(test_batch), function(b) {
                if (!(b %in% names(batch_means))) {
                        return(test_x[test_batch == b, , drop=FALSE]) # unseen batch fallback
                }
                norm_transform(test_x[test_batch == b, , drop=FALSE], b)
        }))
        
        list(train = train_out, test = test_out)
}

impute_train_test <- function(train_x, test_x, method = "mean") {
        
        if (method == "mean"){
                impute_vals <- colMeans(train_x, na.rm = TRUE)
        }else if (method == "minimum"){
                impute_vals <- apply(train_x, 2, function(x) min(x[!is.na(x)]))
        }
        impute_na <- function(x) {
                for (j in seq_len(ncol(x))) {
                        x[is.na(x[, j]), j] <- impute_vals[j]
                }
                x
        }
        
        list(
                train = impute_na(train_x),
                test  = impute_na(test_x)
        )
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

get_r2 <- function(y, yhat){
        r2 <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
        return(r2)
}

batches <- unique(batch)

library(glmnet)

set.seed(666)
hold_out_batch <- sample(batches, 1)


library(xgboost)
cell_types <- unique(deg_list$cell_type)

cv_res_list <- list()
final_mod_res <- list()
glm_list <- list()
plot_list <- list()
for(ct in cell_types){
        print(sprintf("Fitting %s model...", ct))
        results_glm <- c()
        r2_glm_cv <- c()
        results_xgb <- c()
        r2_xgb_cv <- c()
        
        X_cell <- X[batch != hold_out_batch, prot_info$protein[prot_info$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]]
        X_cell_test <- X[batch == hold_out_batch, prot_info$protein[prot_info$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]]

        prot_stats_glm_cv <- data.frame(protein = colnames(X_cell),
                                        times_selected = 0)
        selected_alphas <- numeric(length(batches[batches != hold_out_batch]))
        for (b in batches[batches != hold_out_batch]) {
                test_idx <- batch[batch != hold_out_batch] == b
                train_idx <- !test_idx
                
                X_train <- X_cell[train_idx, ]
                X_test  <- X_cell[test_idx, ]
                
                y_train <- y[batch != hold_out_batch][train_idx]
                y_test  <- y[batch != hold_out_batch][test_idx]
                
                batch_train <- batch[batch != hold_out_batch][train_idx]
                batch_test  <- batch[batch != hold_out_batch][test_idx]
                
                # 1. batch correction (fit on train only)
                bc <- batch_normalize(X_train, batch_train, X_test, batch_test)
                
                # 2. imputation (train only)
                imp <- impute_train_test(bc$train, bc$test, method = "mean")
                
                Xtr <- imp$train
                Xte <- imp$test
                
                ##########################################################################
                # GLM
                ##########################################################################
                
                best_alpha <- get_best_alpha(Xtr, y_train)
                selected_alphas <- c(selected_alphas, best_alpha$best_alpha)
                
                cvfit <- cv.glmnet(
                        x = as.matrix(Xtr),
                        y = y_train,
                        family = "gaussian",
                        alpha = best_alpha$best_alpha,
                        nfolds = 5
                )
                
                pred_glm <- predict(cvfit, newx = as.matrix(Xte), s = "lambda.min")
                
                r2_glm <- get_r2(y_test, pred_glm)
                r2_glm_cv <- c(r2_glm_cv, r2_glm)
                results_glm <- c(results_glm, cor(pred_glm, y_test, use="complete.obs"))
                
                coefs <- coef(cvfit, s = "lambda.min")
                prot_diff_zero <- coefs[, 1][coefs[, 1] != 0]
                prot_diff_zero <- prot_diff_zero[names(prot_diff_zero) != "(Intercept)"]
                
                prot_stats_glm_cv$times_selected[
                        prot_stats_glm_cv$protein %in% names(prot_diff_zero)
                ] <- prot_stats_glm_cv$times_selected[
                        prot_stats_glm_cv$protein %in% names(prot_diff_zero)
                ] + 1
                
                ##########################################################################
                # XGBOOST
                ##########################################################################
                
                dtrain <- xgb.DMatrix(as.matrix(Xtr), label = y_train)
                dtest  <- xgb.DMatrix(as.matrix(Xte), label = y_test)
                
                params <- list(
                        objective = "reg:squarederror",
                        max_depth = 4,
                        eta = 0.05,
                        subsample = 0.8,
                        seed = 123,
                        colsample_bytree = 0.8
                )
                
                xgb_fit <- xgb.train(params, dtrain, nrounds = 200, verbose = 0)
                
                pred_xgb <- predict(xgb_fit, dtest)
                
                r2_xgb <- get_r2(y_test, pred_xgb)
                r2_xgb_cv <- c(r2_xgb_cv, r2_xgb)
                results_xgb <- c(results_xgb, cor(pred_xgb, y_test, use="complete.obs"))
        }
        
        mean(results_glm)
        mean(results_xgb)
        mean(r2_glm_cv)
        mean(r2_xgb_cv)
        
        final_alpha <- as.numeric(
                names(
                        which.max(
                                table(selected_alphas)
                        )
                )
        )
        
        prot_stats_glm_cv$times_selected <- prot_stats_glm_cv$times_selected / length(batches[batches != hold_out_batch])
        
        sel_prots <- prot_stats_glm_cv$protein[prot_stats_glm_cv$times_selected >= .5]
        
        X_cell_bc <- batch_normalize(
                X_cell,
                batch[batch != hold_out_batch],
                X_cell_test,
                batch[batch == hold_out_batch])
        #X_cell_bc <- X_cell_bc$train
        
        # 2. imputation (train only)
        X_cell_imp <- impute_train_test(X_cell_bc$train, X_cell_bc$test, method = "mean")
        
        final_cvfit <- cv.glmnet(
                x = as.matrix(X_cell_imp$train[, sel_prots]),
                y = y[batch != hold_out_batch],
                family = "gaussian",
                alpha = final_alpha,
                nfolds = 5
        )
        
        pred_glm_final <- predict(
                final_cvfit,
                newx = as.matrix(X_cell_imp$test[, sel_prots]),
                s = "lambda.min"
        )
        pred_df <- data.frame(
                age = y[batch == hold_out_batch],
                age_pred = pred_glm_final[, 1]
        )
        
        fit <- lm(age_pred ~ age, data = pred_df)
        
        intercept <- coef(fit)[1]
        slope <- coef(fit)[2]
        r <- cor(pred_df$age, pred_df$age_pred, use = "complete.obs")
        
        eq_label <- sprintf(
                "y = %.2f + %.2fx\nr = %.3f",
                intercept, slope, r
        )
        
        cor_plot <- ggplot(pred_df, aes(x = age, y = age_pred)) +
                geom_point() +
                geom_smooth(method = "lm") +
                ggtitle(ct) +
                annotate(
                        "text",
                        x = Inf, y = -Inf,
                        label = eq_label,
                        hjust = 1.1, vjust = -0.5
                ) +
                theme(
                        title = ggtext::element_markdown(),
                        axis.title.y = ggtext::element_markdown(),
                        panel.background = element_rect(
                                fill = "white",
                                color = "white"
                        ),
                        panel.border = element_rect(
                                colour = "black", fill=NA, linewidth = 1
                        ),
                        panel.grid.major = element_line(colour = "#d4d4d4")
                )
        plot_list[[ct]] <- cor_plot
        glm_cor_final <- cor(y[batch == hold_out_batch], pred_glm_final)
        glm_r2_final <- get_r2(y[batch == hold_out_batch], pred_glm_final)
        
        cv_df <- data.frame(
                cell_type = ct,
                fold =  batches[batches != hold_out_batch],
                glm_cor = results_glm,
                glm_r2 = r2_glm_cv,
                xgb_cor = results_xgb,
                xgb_r2 = r2_xgb_cv
        )
        cv_res_list[[ct]] <- cv_df
        final_mod_perf <- data.frame(cell_type = ct, glm_cor = glm_cor_final, glm_r2 = glm_r2_final)
        final_mod_res[[ct]] <- final_mod_perf
        glm_list[[ct]] <- final_cvfit
}

final_mod_res <- do.call(rbind.data.frame, final_mod_res)

cv_res_list <- do.call(rbind.data.frame, cv_res_list)

cv_res_mean <- cv_res_list %>%
        group_by(cell_type) %>%
        dplyr::summarize(mean_glm_cor = mean(glm_cor),
                         mean_glm_r2 = mean(glm_r2),
                         mean_xgb_cor = mean(xgb_cor),
                         mean_xgb_r2 = mean(xgb_r2))

final_mod_res
cv_res_mean

plot_list

colMeans(is.na(X[, prot_info$protein[prot_info$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == "Microglia"]]]))
colMeans(is.na(X[, prot_info$protein[prot_info$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == "Astrocyte"]]]))

# Train a model on everything to compare disease differentiation without looking 
# at the cell-specific processes

dim(X)

X_train <- X[mdat$Sample_preparation_batch != hold_out_batch, ]
X_test <- X[mdat$Sample_preparation_batch == hold_out_batch, ]


X_all_bc <- batch_normalize(
        X_train,
        mdat$Sample_preparation_batch[
                mdat$Sample_preparation_batch != hold_out_batch
        ],
        X_test,
        mdat$Sample_preparation_batch[
                mdat$Sample_preparation_batch == hold_out_batch
        ]
)

X_all_imp <- impute_train_test(X_all_bc$train, X_all_bc$test, method = "mean")

best_alpha <- get_best_alpha(
        X_all_imp$train,
        mdat$Age[mdat$Sample_preparation_batch != hold_out_batch]
)

allmod_cvfit <- cv.glmnet(
        x = as.matrix(X_all_imp$train),
        y = y[batch != hold_out_batch],
        family = "gaussian",
        alpha = best_alpha$best_alpha,
        nfolds = 5
)

pred_glm_allmod <- predict(
        allmod_cvfit,
        newx = as.matrix(X_all_imp$test),
        s = "lambda.min"
)

plot(mdat$Age[mdat$Sample_preparation_batch == hold_out_batch], pred_glm_allmod[, 1])

print(
        sprintf(
                "Full model R: %s",
                cor(
                        mdat$Age[mdat$Sample_preparation_batch == hold_out_batch],
                        pred_glm_allmod[, 1]
                )
        )
)

print(
        sprintf(
                "Full model R2: %s",
                get_r2(
                        mdat$Age[mdat$Sample_preparation_batch == hold_out_batch],
                        pred_glm_allmod[, 1]
                )
        )
)




# Load complete bader dataset to see what's the outcome in the disease ones
################################################################################

bader_complete <- plotUtils::read_table_fast("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/DA-F08.4_-SEC-pass_v06Sc_ion_LibPGQVal1perc_precdLFQv03s01_prot_preprSc04.tsv", row.names = 1)
rownames(bader_complete) <- bader_complete$protein
bader_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots/bader_mdat_parsed.csv", row.names = 1)

bader_complete <- t(
        bader_complete[,
                match(
                        make.names(bader_mdat$Run_timsTOFMeasurement),
                        make.names(colnames(bader_complete))
                )
        ]
)

bader_complete <- log2(bader_complete)

# Compute acceleration on each cell type and assess in the disease individuals
################################################################################


boxplot(Age ~ Diagnosis_group, data=bader_mdat)

cell_types
accel_plots <- list()
acceleration_df <- list()
disease_cell_assoc <- list()
for (ct in cell_types){
        print(sprintf("Computing acceleration in %s...", ct))
        mod_cell <- glm_list[[ct]]
        mod_genes <- rownames(coef(mod_cell))
        mod_genes <- mod_genes[mod_genes != "(Intercept)"]
        X_cell <- bader_complete[, mod_genes]
        
        samps_train <- rownames(X_cell) %in% rownames(X[mdat$Sample_preparation_batch != hold_out_batch, ])
        X_cell_train <- X_cell[samps_train, ]
        X_cell_test <- X_cell[!samps_train, ]
        
        batch_train <- bader_mdat$Sample_preparation_batch[samps_train]
        batch_test <- bader_mdat$Sample_preparation_batch[!samps_train]
        
        X_cell_bc <- batch_normalize(
                X_cell_train,
                batch_train,
                X_cell_test,
                batch_test
        )
        #X_cell_bc <- X_cell_bc$train
        
        # 2. imputation (train only)
        X_cell_imp <- impute_train_test(
                X_cell_bc$train,
                X_cell_bc$test,
                method = "mean"
        )
        X_cell_imp <- rbind.data.frame(
                X_cell_imp$train,
                X_cell_imp$test
        )
        X_cell_imp <- X_cell_imp[
                match(
                        make.names(bader_mdat$Run_timsTOFMeasurement),
                        rownames(X_cell_imp)
                ),
        ]
        age_pred <- predict(
                mod_cell,
                newx = as.matrix(X_cell_imp),
                s = "lambda.min"
        )
        pred_df <- data.frame(
                sample = bader_mdat$Run_timsTOFMeasurement,
                cell_type = ct,
                age = bader_mdat$Age,
                age_pred = age_pred[, 1],
                diagnosis = bader_mdat$Diagnosis_group,
                batch = bader_mdat$Sample_preparation_batch
        )
        rownames(pred_df) <- 1:nrow(pred_df)
        
        
        ggplot(
                data = pred_df,
                mapping = aes(x = age, y = age_pred, colour = diagnosis)
        ) +
                geom_point()
        pred_df$raw_accel <- pred_df$age_pred - pred_df$age
        
        accel_mod_1 <- lm(
                raw_accel ~ age + batch,
                data = pred_df[pred_df$diagnosis == "Neurological Control", ]
        )
        
        accel_mod_2 <- lm(
                raw_accel ~ splines::ns(age, df=3) + batch,
                data = pred_df[pred_df$diagnosis == "Neurological Control", ]
        )
        
        accel_mod_3 <- lm(
                raw_accel ~ splines::ns(age, df=4) + batch,
                data = pred_df[pred_df$diagnosis == "Neurological Control", ]
        )
        
        accel_mod_4 <- lm(
                raw_accel ~ splines::ns(age, df=5) + batch,
                data = pred_df[pred_df$diagnosis == "Neurological Control", ]
        )
        
        AIC(accel_mod_1, accel_mod_2, accel_mod_3, accel_mod_4)
        
        pred_df$acceleration <- pred_df$raw_accel - predict(accel_mod_1, pred_df)
        accel_plot <- ggplot(
                data = pred_df,
                mapping = aes(x = age, y = acceleration, colour = diagnosis)
        ) +
                geom_point() +
                geom_smooth(method = "lm")
        diseases <- unique(pred_df$diagnosis)
        pred_df$diagnosis <- factor(
                pred_df$diagnosis,
                levels = c("Neurological Control", diseases[diseases != "Neurological Control"])
        )
        #coef(summary(lm(age_pred ~ diagnosis + age + batch, data = pred_df)))
        fit <- lm(
                acceleration ~ diagnosis,
                data = pred_df
        )
        ggplot(pred_df, mapping = aes(x = diagnosis, y = acceleration)) +
                geom_violin()
        
        fit_res <- as.data.frame(coef(summary(fit)))
        fit_res$disease <- rownames(fit_res)
        fit_res$cell_type <- ct
        fit_res <- fit_res[
                fit_res$disease != "(Intercept)",
                c(ncol(fit_res) - 1, ncol(fit_res), 1:(ncol(fit_res) - 2))
        ]
        accel_plots[[ct]] <- accel_plot
        acceleration_df[[ct]] <- pred_df
        disease_cell_assoc[[ct]] <- fit_res
}

acceleration_df <- do.call(rbind.data.frame, acceleration_df)
disease_cell_assoc <- do.call(rbind.data.frame, disease_cell_assoc)
disease_cell_assoc[order(disease_cell_assoc$disease), ]

ctrl_df <- pred_df[pred_df$diagnosis == "Neurological Control", ]
plot(ctrl_df$age, ctrl_df$raw_accel)

plot(density(mdat_all$QAlb_log2, na.rm = TRUE))
acceleration_df$QAlb_log2 <- bader_mdat$QAlb_log2[
        match(acceleration_df$sample, bader_mdat$Run_timsTOFMeasurement)
]

acceleration_df$OCB <- bader_mdat$Oligoclonal_bands[
        match(acceleration_df$sample, bader_mdat$Run_timsTOFMeasurement)
]

dis_pheno_df <- list()
for(ct in cell_types){
        acceleration_df_ct <- acceleration_df[acceleration_df$cell_type == ct, ]
        qalb_cor <- cor(
                acceleration_df_ct$acceleration,
                acceleration_df_ct$QAlb_log2,
                use = "complete.obs"
        )
        qalb_cor_pval <- cor.test(acceleration_df_ct$acceleration,
                 acceleration_df_ct$QAlb_log2,
                 use = "complete.obs")$p.value
        ocb_t_test <- t.test(
                acceleration_df_ct$acceleration[acceleration_df_ct$OCB == "negative"],
                acceleration_df_ct$acceleration[acceleration_df_ct$OCB == "positive"],
                alternative = "less"
        )
        res_df <- data.frame(
                cell_type = ct,
                qalb_cor = qalb_cor,
                qalb_cor_pvalue = qalb_cor_pval,
                ocb_eff_size = ocb_t_test$estimate[2] - ocb_t_test$estimate[1],
                ocb_pval = ocb_t_test$p.value
        )
        dis_pheno_df[[ct]] <- res_df
}
dis_pheno_df <- do.call(rbind.data.frame, dis_pheno_df)


ggplot(
        acceleration_df[
                acceleration_df$diagnosis %in% c("Neurological Control", "Neurodegenerative"),
        ],
        mapping = aes(x = age, y = acceleration, color = diagnosis)
) +
        geom_point() + geom_smooth(method = "lm")

ggplot(
        acceleration_df[
                acceleration_df$diagnosis %in% c("Neurological Control", "Neurodegenerative") & acceleration_df$cell_type == "Astrocyte",
        ],
        mapping = aes(x = age, y = raw_accel, color = diagnosis)
) +
        geom_point() + geom_smooth(method = "lm")


acceleration_celltype_mat <- acceleration_df %>%
        dplyr::select(sample, cell_type, acceleration) %>%
        tidyr::pivot_wider(
                id_cols = sample,
                names_from = cell_type,
                values_from = acceleration,
                names_prefix = "accel_"
        ) %>%
        as.data.frame()

rownames(acceleration_celltype_mat) <- acceleration_celltype_mat$sample

acceleration_celltype_mat <- acceleration_celltype_mat[, 2:ncol(acceleration_celltype_mat)]

accel_pca <- prcomp(acceleration_celltype_mat, scale. = TRUE, center = TRUE)
bader_mdat$sample <- bader_mdat$Run_timsTOFMeasurement



plotUtils::plotPCA(accel_pca, samp_info = bader_mdat, col = "Diagnosis_group")

plotUtils::plotPCA(accel_pca, samp_info = bader_mdat, col = "QAlb_log2")

if(!require("uwot", quietly = TRUE)){
        install.packages("uwot")
}

library(uwot)
library(ggplot2)

umap_res <- uwot::umap(
        accel_pca$x[, 1:4],   # use first PCs (adjust 10–50)
        n_neighbors = 5,
        min_dist = 0.1,
        metric = "euclidean"
)

umap_df <- data.frame(
        sample = rownames(accel_pca$x),
        UMAP1 = umap_res[, 1],
        UMAP2 = umap_res[, 2]
)

umap_df <- dplyr::left_join(
        umap_df,
        bader_mdat,
        by = "sample"
)

diagnosis_df <- as.data.frame.matrix(table(umap_df$sample, umap_df$Diagnosis_group))
colnames(diagnosis_df) <- paste0(make.names(colnames(diagnosis_df)), "_binary")
#diagnosis_df <- data.frame(diagnosis_df)
head(diagnosis_df)
umap_df <- cbind.data.frame(umap_df, diagnosis_df)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Diagnosis_group)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Infectious_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = MS_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Neoplasia_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Neurodegenerative_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Neurological.Control_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Other_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = OtherAutoimmune_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Stroke_binary)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic()


cor(acceleration_celltype_mat)


# Compute acceleration using the complete model
################################################################################

all(make.names(bader_mdat$Run_timsTOFMeasurement) == rownames(bader_complete))
samps_train <- rownames(bader_complete) %in% rownames(X[mdat$Sample_preparation_batch != hold_out_batch, ])

allmod_feats <- rownames(
        coef(allmod_cvfit)
)
allmod_feats <- allmod_feats[allmod_feats != "(Intercept)"]

X_train_4wholemod_pred <- bader_complete[samps_train, allmod_feats]
X_rest_4wholemod_pred <- bader_complete[!samps_train, allmod_feats]


X_bc_4wholemod_pred <- batch_normalize(
        X_train_4wholemod_pred,
        bader_mdat$Sample_preparation_batch[
                match(
                        rownames(X_train_4wholemod_pred),
                        make.names(bader_mdat$Run_timsTOFMeasurement)
                )
        ],
        X_rest_4wholemod_pred,
        bader_mdat$Sample_preparation_batch[
                match(
                        rownames(X_rest_4wholemod_pred),
                        make.names(bader_mdat$Run_timsTOFMeasurement)
                )
        ]
)


X_imp_4wholemod_pred <- impute_train_test(
        X_bc_4wholemod_pred$train,
        X_bc_4wholemod_pred$test,
        method = "mean"
)
X_imp_4wholemod_pred <- rbind.data.frame(
        X_imp_4wholemod_pred$train,
        X_imp_4wholemod_pred$test
)

X_imp_4wholemod_pred <- X_imp_4wholemod_pred[
        match(
                make.names(bader_mdat$Run_timsTOFMeasurement),
                rownames(X_imp_4wholemod_pred)
        ),
]
age_pred_wholemod_pred <- predict(
        allmod_cvfit,
        newx = as.matrix(X_imp_4wholemod_pred),
        s = "lambda.min"
)

pred_df_wholemod_pred <- data.frame(
        sample = bader_mdat$Run_timsTOFMeasurement,
        age = bader_mdat$Age,
        age_pred = age_pred_wholemod_pred[, 1],
        diagnosis = bader_mdat$Diagnosis_group,
        batch = bader_mdat$Sample_preparation_batch
)

pred_df_wholemod_pred$raw_accel <- pred_df_wholemod_pred$age_pred - pred_df_wholemod_pred$age


accel_mod_wholemod_pred <- lm(
        raw_accel ~ age + batch,
        data = pred_df_wholemod_pred[pred_df_wholemod_pred$diagnosis == "Neurological Control", ]
)

pred_df_wholemod_pred$acceleration <- pred_df_wholemod_pred$raw_accel - predict(accel_mod_wholemod_pred, pred_df_wholemod_pred)

pred_df_wholemod_pred$diagnosis <- factor(
        pred_df_wholemod_pred$diagnosis,
        levels = c(
                "Neurological Control",
                unique(pred_df_wholemod_pred$diagnosis)[unique(pred_df_wholemod_pred$diagnosis) != "Neurological Control"]
        )
)

summary(lm(acceleration ~ diagnosis, data = pred_df_wholemod_pred))


# See if cell-specific acceleration association to disease remains after
# incuding the global acceleration as covariate
acceleration_df$acceleration_global_model <- pred_df_wholemod_pred$acceleration[
        match(acceleration_df$sample, pred_df_wholemod_pred$sample)
]

disease_cell_assoc_wholemod <- list()
for (ct in cell_types){
        acceleration_df_ct <- acceleration_df[acceleration_df$cell_type == ct, ]
        #acceleration_df_ct$diagnosis <- factor(
        #        pred_df$diagnosis,
        #        levels = c("Neurological Control", diseases[diseases != "Neurological Control"])
        #)
        #coef(summary(lm(age_pred ~ diagnosis + age + batch, data = pred_df)))
        fit <- lm(
                acceleration ~ diagnosis + acceleration_global_model + QAlb_log2,
                data = acceleration_df_ct
        )
        
        fit_res <- as.data.frame(coef(summary(fit)))
        fit_res$disease <- rownames(fit_res)
        fit_res$cell_type <- ct
        fit_res <- fit_res[
                fit_res$disease != "(Intercept)",
                c(ncol(fit_res) - 1, ncol(fit_res), 1:(ncol(fit_res) - 2))
        ]
        disease_cell_assoc_wholemod[[ct]] <- fit_res
}

disease_cell_assoc_wholemod <- do.call(
        rbind.data.frame, disease_cell_assoc_wholemod
)

disease_cell_assoc_wholemod[order(disease_cell_assoc_wholemod$disease), ]


# Repeat, but with MS disease subtypes

acceleration_df$diagnosis_detailed <- as.character(acceleration_df$diagnosis)

acceleration_df$diagnosis_detailed[
        acceleration_df$diagnosis_detailed == "MS"
] <- bader_mdat$RMS_PMS[bader_mdat$Diagnosis_group == "MS"]

#acceleration_df$diagnosis_detailed[
#        acceleration_df$diagnosis_detailed == "OtherAutoimmune"
#] <- bader_mdat$Diagnosis_group_autoimmune_split_short[bader_mdat$Diagnosis_group == "OtherAutoimmune"]


table(bader_mdat$RMS_PMS, bader_mdat$Diagnosis_group)
table(acceleration_df$diagnosis_detailed)


acceleration_df$diagnosis_detailed <- factor(
        acceleration_df$diagnosis_detailed,
        levels = c(
                "Neurological Control",
                unique(acceleration_df$diagnosis_detailed)[
                        unique(acceleration_df$diagnosis_detailed) != "Neurological Control"
                ]
        )
)

disease_cell_assoc_wholemod_subdiag <- list()
for (ct in cell_types){
        acceleration_df_ct <- acceleration_df[acceleration_df$cell_type == ct, ]
        #acceleration_df_ct$diagnosis <- factor(
        #        pred_df$diagnosis,
        #        levels = c("Neurological Control", diseases[diseases != "Neurological Control"])
        #)
        #coef(summary(lm(age_pred ~ diagnosis + age + batch, data = pred_df)))
        fit <- lm(
                acceleration ~ diagnosis_detailed + acceleration_global_model + QAlb_log2,
                data = acceleration_df_ct
        )
        
        fit_res <- as.data.frame(coef(summary(fit)))
        fit_res$disease <- rownames(fit_res)
        fit_res$cell_type <- ct
        fit_res <- fit_res[
                fit_res$disease != "(Intercept)",
                c(ncol(fit_res) - 1, ncol(fit_res), 1:(ncol(fit_res) - 2))
        ]
        disease_cell_assoc_wholemod_subdiag[[ct]] <- fit_res
}

disease_cell_assoc_wholemod_subdiag <- do.call(
        rbind.data.frame, disease_cell_assoc_wholemod_subdiag
)

disease_cell_assoc_wholemod_subdiag[order(disease_cell_assoc_wholemod_subdiag$disease), ]

ggplot(acceleration_df, mapping = aes(x = diagnosis_detailed, y = QAlb_log2)) +
        geom_boxplot()


disease_cell_assoc_subdiag <- list()
for (ct in cell_types){
        acceleration_df_ct <- acceleration_df[acceleration_df$cell_type == ct, ]
        #acceleration_df_ct$diagnosis <- factor(
        #        pred_df$diagnosis,
        #        levels = c("Neurological Control", diseases[diseases != "Neurological Control"])
        #)
        #coef(summary(lm(age_pred ~ diagnosis + age + batch, data = pred_df)))
        fit <- lm(
                acceleration ~ diagnosis_detailed,
                data = acceleration_df_ct
        )
        
        fit_res <- as.data.frame(coef(summary(fit)))
        fit_res$disease <- rownames(fit_res)
        fit_res$cell_type <- ct
        fit_res <- fit_res[
                fit_res$disease != "(Intercept)",
                c(ncol(fit_res) - 1, ncol(fit_res), 1:(ncol(fit_res) - 2))
        ]
        disease_cell_assoc_subdiag[[ct]] <- fit_res
}

disease_cell_assoc_subdiag <- do.call(
        rbind.data.frame, disease_cell_assoc_subdiag
)

# Do some plotting of the results
################################################################################

ct_order <- cv_res_mean$cell_type[order(cv_res_mean$mean_glm_cor)]
cv_res_mean$cell_type <- factor(cv_res_mean$cell_type, levels = ct_order)
cv_res_list$cell_type <- factor(cv_res_list$cell_type, levels = ct_order)

p_perf <- ggplot(cv_res_mean, aes(x = cell_type, y = mean_glm_cor)) +
        geom_col(fill = "#3B6E8C", width = 0.6) +
        geom_jitter(data = cv_res_list, aes(y = glm_cor), width = 0.1,
                    alpha = 0.5, size = 1.8, inherit.aes = TRUE) +
        coord_flip() +
        labs(x = NULL, y = expression("Cross-validated " * R),
             title = "Clock performance by cell type",
             subtitle = "7-fold batch-blocked CV; points = individual folds") +
        theme_classic(base_size = 13)

p_perf

ggsave(filename = file.path(out_dir, "cv_performance.pdf"), plot = p_perf,
       height = 3, width = 4)

if (!require("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

p_grid <- wrap_plots(plot_list, ncol = 4) +
        plot_annotation(title = "Predicted vs. chronological age (held-out batch)")

p_grid

ggsave(filename = file.path(out_dir, "heldout_performance.pdf"), plot = p_grid,
       height = 5, width = 10)

perf_df <- cv_res_mean %>%
        select(cell_type, CV = mean_glm_cor) %>%
        left_join(
                final_mod_res %>% select(cell_type, HeldOut = lambda.min),
                by = "cell_type"
        ) %>%
        pivot_longer(
                cols = c(CV, HeldOut),
                names_to = "evaluation",
                values_to = "R"
        )

cv_and_ho_plot <- ggplot() +
        geom_col(
                data = perf_df,
                aes(cell_type, R, fill = evaluation),
                position = position_dodge(0.7),
                width = 0.6
        ) +
        geom_jitter(
                data = cv_res_list,
                aes(cell_type, glm_cor),
                width = 0.08,
                alpha = 0.5,
                size = 1.8
        ) +
        coord_flip() +
        theme_classic()

ggsave(file = file.path(out_dir, "cv_and_ho_plot.pdf"), plot = cv_and_ho_plot, height = 3, width = 4)

# Heatmap disease x cell effect size
heatmap_df <- disease_cell_assoc_wholemod %>%
        filter(!disease %in% c("acceleration_global_model", "QAlb_log2")) %>%
        mutate(
                disease_clean = gsub("^diagnosis", "", disease),
                p_adj = p.adjust(`Pr(>|t|)`, method = "BH"),
                sig_label = case_when(
                        p_adj < 0.001 ~ "***",
                        p_adj < 0.01  ~ "**",
                        p_adj < 0.05  ~ "*",
                        TRUE ~ ""
                )
        )

# order axes by overall effect magnitude for a cleaner layout
disease_order <- heatmap_df %>%
        group_by(disease_clean) %>%
        summarize(m = mean(abs(Estimate))) %>%
        arrange(m) %>% pull(disease_clean)

celltype_order <- heatmap_df %>%
        group_by(cell_type) %>%
        summarize(m = mean(abs(Estimate))) %>%
        arrange(m) %>% pull(cell_type)

heatmap_df <- heatmap_df %>%
        mutate(
                disease_clean = factor(disease_clean, levels = disease_order),
                cell_type = factor(cell_type, levels = celltype_order)
        )

p_heatmap <- ggplot(heatmap_df, aes(x = disease_clean, y = cell_type, fill = Estimate)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_label), vjust = 0.75, size = 4) +
        scale_fill_gradient2(
                low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                name = "Acceleration\n(years)"
        ) +
        labs(x = NULL, y = NULL#,
             #title = "Cell-type-specific aging acceleration by diagnosis",
             #subtitle = "Adjusted for global acceleration & QAlb; *p_adj<0.05, **<0.01, ***<0.001"
        ) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_heatmap

# Same, but with MS subtypes:

# Heatmap disease x cell effect size
heatmap_df <- disease_cell_assoc_wholemod_subdiag %>%
        filter(!disease %in% c("acceleration_global_model", "QAlb_log2")) %>%
        mutate(disease_clean = gsub("^diagnosis_detailed", "", disease)) %>%
        group_by(cell_type) %>%
        mutate(
                p_adj = p.adjust(`Pr(>|t|)`, method = "BH"),
                sig_label = case_when(
                        p_adj < 0.001 ~ "***",
                        p_adj < 0.01  ~ "**",
                        p_adj < 0.05  ~ "*",
                        TRUE ~ ""
                )
        ) %>%
        ungroup()

heatmap_df$disease_clean <- gsub("MS_", "", heatmap_df$disease_clean)

# order axes by overall effect magnitude for a cleaner layout
disease_order <- heatmap_df %>%
        group_by(disease_clean) %>%
        summarize(m = mean(abs(Estimate))) %>%
        arrange(m) %>% pull(disease_clean)

celltype_order <- heatmap_df %>%
        group_by(cell_type) %>%
        summarize(m = mean(abs(Estimate))) %>%
        arrange(m) %>% pull(cell_type)

heatmap_df <- heatmap_df %>%
        mutate(
                disease_clean = factor(disease_clean, levels = disease_order),
                cell_type = factor(cell_type, levels = celltype_order)
        )

p_heatmap <- ggplot(heatmap_df, aes(x = disease_clean, y = cell_type, fill = Estimate)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_label), vjust = 0.75, size = 4) +
        scale_fill_gradient2(
                low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                name = "Acceleration\n(years)"
        ) +
        labs(x = NULL, y = NULL#,
             #title = "Cell-type-specific aging acceleration by diagnosis",
             #subtitle = "Adjusted for global acceleration & QAlb; *p_adj<0.05, **<0.01, ***<0.001"
        ) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_heatmap

ggsave(filename = file.path(out_dir, "cell_disease_hm.pdf"), plot = p_heatmap, height = 3, width = 5)


# Adjusted vs unadjusted significance:

if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

model1_df <- disease_cell_assoc_subdiag %>%
        filter(!disease %in% c("acceleration_global_model", "QAlb_log2")) %>%
        select(disease, cell_type, p_model1 = `Pr(>|t|)`, est_model1 = Estimate, se_model1 = `Std. Error`)

model1_df$disease_clean <- gsub("diagnosis|diagnosis_detailed", "", model1_df$disease)

compare_df <- heatmap_df %>%
        select(disease, cell_type, p_model3 = `Pr(>|t|)`, est_model3 = Estimate, se_model3 = `Std. Error`, disease_clean) %>%
        inner_join(model1_df, by = c("disease_clean", "cell_type")) %>%
        mutate(
                neglog10_p1 = -log10(p_model1),
                neglog10_p3 = -log10(p_model3),
                label = paste0(cell_type, "/", disease_clean)
        )

p_compare <- ggplot(compare_df, aes(x = neglog10_p1, y = neglog10_p3, color = cell_type)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey70") +
        geom_point(size = 3, alpha = 0.85) +
        geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 25, show.legend = FALSE) +
        labs(
                x = expression(-log[10]*"(p) unadjusted"),
                y = expression(-log[10]*"(p) adjusted for global accel. + QAlb"),
                title = "Robustness of cell-type-specific disease associations",
                color = "Cell type"
        ) +
        theme_classic(base_size = 13)

p_compare


p_compare <- ggplot(compare_df, aes(x = est_model1, y = est_model3, color = cell_type)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
        geom_point(size = 3, alpha = 0.85) +
        geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 25, show.legend = FALSE) +
        labs(
                x = "Acceleration effect (years), unadjusted",
                y = "Acceleration effect (years), adjusted for global accel. + QAlb",
                title = "Cell-type acceleration effects are preserved after adjustment",
                color = "Cell type"
        ) +
        theme_classic(base_size = 13)

p_compare

# RRMS vs PMS
ms_subtype_df <- disease_cell_assoc_wholemod_subdiag %>%
        filter(disease %in% c("diagnosis_detailedRRMS", "diagnosis_detailedPMS")) %>%
        mutate(
                subtype = factor(gsub("diagnosis_detailed", "", disease), levels = c("RRMS", "PMS")),
                ci_lower = Estimate - 1.96 * `Std. Error`,
                ci_upper = Estimate + 1.96 * `Std. Error`,
                p_adj = p.adjust(`Pr(>|t|)`, method = "BH")
        )

p_subtype <- ggplot(ms_subtype_df, aes(x = cell_type, y = Estimate, color = subtype)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper),
                        position = position_dodge(width = 0.5), size = 0.9) +
        coord_flip() +
        labs(x = NULL, y = "Acceleration (years)", color = "MS subtype",
             #title = "Cell-type acceleration in MS: RRMS vs. progressive MS",
             #subtitle = "Adjusted for global acceleration & QAlb; error bars = 95% CI"
             ) +
        theme_classic(base_size = 13)

p_subtype

ggsave(plot = p_subtype, filename = file.path(out_dir, "PMS_vs_RRMS.pdf"), height = 4, width = 4)



















diseases <- unique(acceleration_df$diagnosis)
diseases <- diseases[diseases != "Neurological Control"]
for(d in diseases){
        acceleration_df_d <- acceleration_df[
                acceleration_df$diagnosis %in% c("Neurological Control", d),
        ]
        acceleration_df_d$diagnosis <- factor(
                acceleration_df_d$diagnosis,
                levels = c("Neurological Control", d)
        )
        fit <- lmerTest::lmer(
                acceleration ~ cell_type + (1 | sample),
                data = acceleration_df_d
        )
        coef(summary(fit))
}











library(lme4)

lme4::lmer()

cell_types <- unique(deg_list$cell_type)

cv_res_list <- list()
final_mod_res <- list()
glm_list <- list()
for(ct in cell_types){
        print(sprintf("Fitting %s model...", ct))
        results_glm <- c()
        r2_glm_cv <- c()
        results_xgb <- c()
        r2_xgb_cv <- c()
        
        X_cell <- X[, prot_info$protein[prot_info$ensembl_gene_id %in% deg_list$gene[deg_list$cell_type == ct]]]
        
        prot_stats_glm_cv <- data.frame(protein = colnames(X_cell),
                                        times_selected = 0)
        selected_alphas <- numeric(length(batches))
        for (b in batches) {
                
                test_idx <- batch == b
                train_idx <- !test_idx
                
                X_train <- X_cell[train_idx, ]
                X_test  <- X_cell[test_idx, ]
                
                y_train <- y[train_idx]
                y_test  <- y[test_idx]
                
                batch_train <- batch[train_idx]
                batch_test  <- batch[test_idx]
                
                # 1. batch correction (fit on train only)
                bc <- batch_normalize(X_train, batch_train, X_test, batch_test)
                
                # 2. imputation (train only)
                imp <- impute_train_test(bc$train, bc$test, method = "mean")
                
                Xtr <- imp$train
                Xte <- imp$test
                
                ##########################################################################
                # GLM
                ##########################################################################
                
                best_alpha <- get_best_alpha(Xtr, y_train)
                selected_alphas <- c(selected_alphas, best_alpha$best_alpha)
                
                cvfit <- cv.glmnet(
                        x = as.matrix(Xtr),
                        y = y_train,
                        family = "gaussian",
                        alpha = best_alpha$best_alpha,
                        nfolds = 5
                )
                
                pred_glm <- predict(cvfit, newx = as.matrix(Xte), s = "lambda.min")
                
                r2_glm <- get_r2(y_test, pred_glm)
                r2_glm_cv <- c(r2_glm_cv, r2_glm)
                results_glm <- c(results_glm, cor(pred_glm, y_test, use="complete.obs"))
                
                coefs <- coef(cvfit, s = "lambda.min")
                prot_diff_zero <- coefs[, 1][coefs[, 1] != 0]
                prot_diff_zero <- prot_diff_zero[names(prot_diff_zero) != "(Intercept)"]
                
                prot_stats_glm_cv$times_selected[
                        prot_stats_glm_cv$protein %in% names(prot_diff_zero)
                ] <- prot_stats_glm_cv$times_selected[
                        prot_stats_glm_cv$protein %in% names(prot_diff_zero)
                ] + 1
                
                ##########################################################################
                # XGBOOST
                ##########################################################################
                
                dtrain <- xgb.DMatrix(as.matrix(Xtr), label = y_train)
                dtest  <- xgb.DMatrix(as.matrix(Xte), label = y_test)
                
                params <- list(
                        objective = "reg:squarederror",
                        max_depth = 4,
                        eta = 0.05,
                        subsample = 0.8,
                        seed = 123,
                        colsample_bytree = 0.8
                )
                
                xgb_fit <- xgb.train(params, dtrain, nrounds = 200, verbose = 0)
                
                pred_xgb <- predict(xgb_fit, dtest)
                
                r2_xgb <- get_r2(y_test, pred_xgb)
                r2_xgb_cv <- c(r2_xgb_cv, r2_xgb)
                results_xgb <- c(results_xgb, cor(pred_xgb, y_test, use="complete.obs"))
        }
        
        mean(results_glm)
        mean(results_xgb)
        mean(r2_glm_cv)
        mean(r2_xgb_cv)
        
        final_alpha <- as.numeric(
                names(
                        which.max(
                                table(selected_alphas)
                        )
                )
        )
        
        prot_stats_glm_cv$times_selected <- prot_stats_glm_cv$times_selected / length(batches)
        
        sel_prots <- prot_stats_glm_cv$protein[prot_stats_glm_cv$times_selected >= .5]
        
        X_cell_bc <- batch_normalize(
                X_cell,
                batch,
                X_cell,
                batch)
        #X_cell_bc <- X_cell_bc$train
        
        # 2. imputation (train only)
        X_cell_imp <- impute_train_test(X_cell_bc$train, X_cell_bc$train, method = "mean")
        
        final_cvfit <- cv.glmnet(
                x = as.matrix(X_cell_imp$train[, sel_prots]),
                y = y,
                family = "gaussian",
                alpha = final_alpha,
                nfolds = 5
        )
        
        pred_glm_final <- predict(
                final_cvfit,
                newx = as.matrix(X_cell_imp$train[, sel_prots]),
                s = "lambda.min"
        )
        
        glm_cor_final <- cor(y, pred_glm_final)
        glm_r2_final <- get_r2(y, pred_glm_final)
        
        cv_df <- data.frame(
                cell_type = ct,
                fold =  batches,
                glm_cor = results_glm,
                glm_r2 = r2_glm_cv,
                xgb_cor = results_xgb,
                xgb_r2 = r2_xgb_cv
        )
        cv_res_list[[ct]] <- cv_df
        final_mod_perf <- data.frame(cell_type = ct, glm_cor = glm_cor_final, glm_r2 = glm_r2_final)
        final_mod_res[[ct]] <- final_mod_perf
        glm_list[[ct]] <- final_cvfit
}

final_mod_res <- do.call(rbind.data.frame, final_mod_res)

cv_res_list <- do.call(rbind.data.frame, cv_res_list)

cv_res_mean <- cv_res_list %>%
        group_by(cell_type) %>%
        dplyr::summarize(mean_glm_cor = mean(glm_cor),
                         mean_glm_r2 = mean(glm_r2),
                         mean_xgb_cor = mean(xgb_cor),
                         mean_xgb_r2 = mean(xgb_r2))

final_mod_res
cv_res_mean

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

# Remove the effect of blood-CSF barrier permeability (QAlb) from every protein
# before fitting/applying the clocks (see qalb_residualize)
residualize_qalb <- TRUE

mdat_runs <- mdat_all[!is.na(mdat_all$Run_timsTOFMeasurement), ]
mdat_runs <- mdat_runs[!duplicated(mdat_runs$Run_timsTOFMeasurement), ]
qalb_lookup <- setNames(mdat_runs$QAlb_log2, make.names(mdat_runs$Run_timsTOFMeasurement))
age_lookup <- setNames(mdat_runs$Age, make.names(mdat_runs$Run_timsTOFMeasurement))


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
        
        # Return rows in the input order (the rbind above groups them by batch,
        # which misaligns them with the outcome vector when samples are not
        # sorted by batch)
        list(train = train_out[rownames(train_x), , drop = FALSE],
             test = test_out[rownames(test_x), , drop = FALSE])
}

# Residualize each protein on QAlb_log2. Slopes are estimated in the training
# samples with age as covariate (so the age signal correlated with QAlb is not
# removed), and only the QAlb term is subtracted, in train and test, so
# prediction needs QAlb but not age. Samples are matched by rownames to
# qalb_lookup/age_lookup; missing QAlb is set to the training median.
qalb_residualize <- function(train_x, test_x){
        if (!residualize_qalb){
                return(list(train = train_x, test = test_x))
        }
        q_tr <- unname(qalb_lookup[rownames(train_x)])
        a_tr <- unname(age_lookup[rownames(train_x)])
        q_te <- unname(qalb_lookup[rownames(test_x)])
        ok_samp <- !is.na(q_tr) & !is.na(a_tr)
        q_med <- median(q_tr, na.rm = TRUE)
        q_mu <- mean(q_tr[ok_samp])
        betas <- apply(train_x, 2, function(p){
                ok <- ok_samp & !is.na(p)
                if (sum(ok) < 10) return(0)
                lm.fit(cbind(1, q_tr[ok], a_tr[ok]), p[ok])$coefficients[2]
        })
        betas[is.na(betas)] <- 0
        q_tr[is.na(q_tr)] <- q_med
        q_te[is.na(q_te)] <- q_med
        list(train = train_x - outer(q_tr - q_mu, betas),
             test = test_x - outer(q_te - q_mu, betas),
             betas = betas)
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
alpha_list <- list()
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
                bc <- qalb_residualize(bc$train, bc$test)
                
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
        X_cell_bc <- qalb_residualize(X_cell_bc$train, X_cell_bc$test)
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
        alpha_list[[ct]] <- final_alpha
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

X_all_bc <- qalb_residualize(X_all_bc$train, X_all_bc$test)

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
allmod_alpha <- best_alpha$best_alpha

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
        X_cell_bc <- qalb_residualize(X_cell_bc$train, X_cell_bc$test)
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


X_bc_4wholemod_pred <- qalb_residualize(
        X_bc_4wholemod_pred$train,
        X_bc_4wholemod_pred$test
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


################################################################################
################################################################################
# ADDITIONAL ANALYSES FOR THE WashU (Cruchaga lab) MEETING                     #
################################################################################
################################################################################
# Sections:
#   A. Out-of-batch (leakage-free) clock predictions and accelerations
#   B. Aging-axis alignment: does disease push each cell-type proteome along
#      its own aging direction?
#   C. Pseudo-temporal dynamics from cross-sectional data
#      (diagnosis x age divergence, abnormality-frequency ordering,
#      clinical/progression variables)
#   D. Age-gappers: disease enrichment in extreme agers and multi-cell aging
#   E. Within age-gappers: do the clock proteins separate disease from
#      healthy agers? (cross-sectional proof of concept of the UKBB analysis)
#   F. Added value of the cell-type acceleration profile for diagnosis
#   G. Sex differences in cell-type acceleration
#   H. Transferability to the Knight ADRC SomaScan 7k platform
#   I. Longitudinal samples in the Bader dataset (if any)
################################################################################

if (!require("lme4", quietly = TRUE)) install.packages("lme4")
library(ggplot2)

meet_dir <- file.path(out_dir, "meeting_analyses")
plotUtils::create_dir_if_not(meet_dir)

min_n_diag <- 15        # minimum n per diagnosis group to be modelled
z_thr_gap <- 1.5        # |z| threshold defining extreme agers (age-gappers)
z_thr_gap_soft <- 1     # softer threshold, used where n would be too small

somascan7k_f <- "/Users/guillem.santamaria/Downloads/csf_soma7k_pQTL_aptamer_info.xlsx"
bader_full_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/DA-F08.4_-SEC-pass_v06Sc_ion_LibPGQVal1perc_precdLFQv03s01_prot_preprSc04.tsv"

# Functions
################################################################################

get_clock_prots <- function(mod){
        p <- rownames(coef(mod))
        p[p != "(Intercept)"]
}

# Re-order the output of batch_normalize (which groups rows by batch) back to
# the original row order
batch_normalize_keep_order <- function(train_x, train_batch, test_x, test_batch){
        bc <- batch_normalize(train_x, train_batch, test_x, test_batch)
        list(train = bc$train[rownames(train_x), , drop = FALSE],
             test = bc$test[rownames(test_x), , drop = FALSE])
}

# Leave-one-batch-out predictions: every sample (control or patient) is
# predicted by a model that never saw any sample from its batch.
# X_bc: batch-normalized matrix (samples x proteins) of all samples.
lobo_clock_predict <- function(prots, alpha, X_bc, age, batch_vec, is_train,
                               nfolds = 5){
        prots <- intersect(prots, colnames(X_bc))
        pred <- setNames(rep(NA_real_, nrow(X_bc)), rownames(X_bc))
        for (b in unique(batch_vec)){
                tr <- is_train & batch_vec != b & !is.na(age)
                te <- batch_vec == b
                keep <- colSums(!is.na(X_bc[tr, prots, drop = FALSE])) > 2
                res <- qalb_residualize(
                        X_bc[tr, prots[keep], drop = FALSE],
                        X_bc[te, prots[keep], drop = FALSE]
                )
                imp <- impute_train_test(res$train, res$test, method = "mean")
                fit <- cv.glmnet(
                        x = as.matrix(imp$train),
                        y = age[tr],
                        family = "gaussian",
                        alpha = alpha,
                        nfolds = nfolds
                )
                pred[te] <- predict(
                        fit, newx = as.matrix(imp$test), s = "lambda.min"
                )[, 1]
        }
        pred
}

# Age-bias correction fitted on controls only. Falls back to an age-only
# correction if some batch has no controls.
age_bias_correct <- function(raw_accel, age, batch_vec, is_ctrl){
        df <- data.frame(raw_accel = raw_accel, age = age, batch = batch_vec)
        if (all(unique(batch_vec) %in% unique(batch_vec[is_ctrl]))){
                fit <- lm(raw_accel ~ age + batch, data = df[is_ctrl, ])
        }else{
                fit <- lm(raw_accel ~ age, data = df[is_ctrl, ])
        }
        df$raw_accel - predict(fit, df)
}

# Fit the same linear model on every clock and return the tidy coefficients
assoc_by_clock <- function(df, rhs, clocks = unique(df$cell_type),
                           drop_terms = "(Intercept)"){
        res <- list()
        for (ct in clocks){
                fit <- lm(
                        as.formula(paste("acceleration ~", rhs)),
                        data = df[df$cell_type == ct, ]
                )
                cf <- as.data.frame(coef(summary(fit)))
                colnames(cf) <- c("estimate", "se", "t", "p")
                cf$term <- rownames(cf)
                cf$cell_type <- ct
                res[[ct]] <- cf[!cf$term %in% drop_terms, ]
        }
        res <- do.call(rbind.data.frame, res)
        rownames(res) <- NULL
        res
}

sig_stars <- function(p){
        dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05 ~ "*", TRUE ~ "")
}

# Rank-based (Mann-Whitney) AUC
auc_mw <- function(y, score){
        ok <- !is.na(score) & !is.na(y)
        y <- as.integer(y[ok])
        r <- rank(score[ok])
        n1 <- sum(y == 1)
        n0 <- sum(y == 0)
        (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_strat_folds <- function(y, k){
        f <- integer(length(y))
        for (lv in unique(y)){
                idx <- which(y == lv)
                f[idx] <- sample(rep_len(seq_len(k), length(idx)))
        }
        f
}

# Repeated stratified CV AUC of a binary outcome.
#   X:      protein matrix (or NULL), fitted with penalized logistic regression
#   covars: data.frame of covariates, left unpenalized (or fitted with glm when
#           X is NULL)
cv_auc <- function(y, X = NULL, covars = NULL, k = 5, n_rep = 5, alpha = 0.5){
        if (!is.null(covars)){
                covars <- droplevels(covars)
                covars <- covars[, sapply(covars, function(v) length(unique(v)) > 1),
                                 drop = FALSE]
        }
        aucs <- c()
        for (r in seq_len(n_rep)){
                folds <- make_strat_folds(y, k)
                score <- rep(NA_real_, length(y))
                for (f in seq_len(k)){
                        tr <- folds != f
                        te <- folds == f
                        if (is.null(X)){
                                d <- data.frame(y = y, covars)
                                fit <- glm(y ~ ., data = d[tr, ], family = binomial)
                                score[te] <- predict(fit, d[te, ], type = "link")
                        }else{
                                keep <- colSums(!is.na(X[tr, , drop = FALSE])) > 0
                                imp <- impute_train_test(
                                        X[tr, keep, drop = FALSE],
                                        X[te, keep, drop = FALSE]
                                )
                                Xtr <- as.matrix(imp$train)
                                Xte <- as.matrix(imp$test)
                                pf <- rep(1, ncol(Xtr))
                                if (!is.null(covars)){
                                        C <- model.matrix(~ ., covars)[, -1, drop = FALSE]
                                        Xtr <- cbind(Xtr, C[tr, , drop = FALSE])
                                        Xte <- cbind(Xte, C[te, , drop = FALSE])
                                        pf <- c(pf, rep(0, ncol(C)))
                                }
                                fit <- cv.glmnet(
                                        Xtr, y[tr], family = "binomial",
                                        alpha = alpha, penalty.factor = pf,
                                        nfolds = min(5, min(table(y[tr])))
                                )
                                score[te] <- predict(fit, newx = Xte, s = "lambda.min")[, 1]
                        }
                }
                aucs <- c(aucs, auc_mw(y, score))
        }
        c(auc = mean(aucs), auc_sd = sd(aucs))
}

# Sample-level metadata and batch-normalized data of the whole cohort
################################################################################
samp_df <- data.frame(
        sample = bader_mdat$Run_timsTOFMeasurement,
        sample_mn = make.names(bader_mdat$Run_timsTOFMeasurement),
        age = bader_mdat$Age,
        sex = factor(bader_mdat$Sex),
        diagnosis = bader_mdat$Diagnosis_group,
        batch = as.character(bader_mdat$Sample_preparation_batch),
        QAlb_log2 = bader_mdat$QAlb_log2
)
samp_df$diagnosis_detailed <- ifelse(
        samp_df$diagnosis == "MS" & !is.na(bader_mdat$RMS_PMS),
        bader_mdat$RMS_PMS,
        samp_df$diagnosis
)
samp_df$is_ctrl <- samp_df$diagnosis == "Neurological Control"
ctrl_lab <- "Neurological Control"
diag_det_levels <- c(
        ctrl_lab,
        sort(setdiff(unique(samp_df$diagnosis_detailed), ctrl_lab))
)
samp_df$diagnosis_detailed <- factor(samp_df$diagnosis_detailed, levels = diag_det_levels)
samp_df$diagnosis <- factor(
        samp_df$diagnosis,
        levels = c(ctrl_lab, sort(setdiff(unique(samp_df$diagnosis), ctrl_lab)))
)

stopifnot(all(rownames(bader_complete) == samp_df$sample_mn))

# Diagnoses with enough samples to model
diag_tab <- table(samp_df$diagnosis_detailed)
diag_model <- setdiff(names(diag_tab)[diag_tab >= min_n_diag], ctrl_lab)
diag_tab

# Samples used for training the clocks (controls)
is_train_ctrl <- rownames(bader_complete) %in% rownames(X)

# Batch normalization of ALL proteins with location/scale estimated in the
# controls of each batch. This is label-free (ages are not used), so all
# batches can contribute without leaking the outcome.
bader_bc <- batch_normalize_keep_order(
        bader_complete[samp_df$is_ctrl, ],
        samp_df$batch[samp_df$is_ctrl],
        bader_complete,
        samp_df$batch
)$test

# QAlb-residualized version (slopes from the training controls), used in the
# protein-level analyses (B, E). Clock predictions residualize within each
# training fold instead (lobo_clock_predict).
bader_bc_q <- qalb_residualize(bader_bc[is_train_ctrl, ], bader_bc)$test

################################################################################
# A. Out-of-batch clock predictions                                           #
################################################################################
# In the acceleration analysis above, controls of the training batches are
# predicted by models trained on themselves (in-sample), while patients and the
# held-out batch are out-of-sample. In-sample predictions are less regressed to
# the mean, so the age-bias correction fitted in controls does not transfer
# exactly to patients, and diagnoses with different age distributions
# (e.g. MS younger, neurodegenerative older) can get spurious accelerations.
# Here every sample is predicted by a model that never saw its batch, keeping
# the protein sets and alphas of the final clocks.

clock_feats <- c(
        lapply(glm_list, get_clock_prots),
        list(Global = intersect(colnames(X), colnames(bader_bc)))
)
clock_alphas <- c(alpha_list, list(Global = allmod_alpha))

set.seed(123)
oof_pred <- sapply(
        names(clock_feats),
        function(ct){
                print(sprintf("Out-of-batch predictions: %s...", ct))
                lobo_clock_predict(
                        clock_feats[[ct]],
                        clock_alphas[[ct]],
                        bader_bc,
                        samp_df$age,
                        samp_df$batch,
                        is_train_ctrl
                )
        }
)

acceleration_oof_df <- list()
for (ct in colnames(oof_pred)){
        df_ct <- samp_df
        df_ct$cell_type <- ct
        df_ct$age_pred <- oof_pred[, ct]
        df_ct$raw_accel <- df_ct$age_pred - df_ct$age
        df_ct$acceleration <- age_bias_correct(
                df_ct$raw_accel, df_ct$age, df_ct$batch, df_ct$is_ctrl
        )
        # Express as SDs of the control distribution, to compare clocks
        df_ct$accel_z <- df_ct$acceleration / sd(df_ct$acceleration[df_ct$is_ctrl], na.rm = TRUE)
        acceleration_oof_df[[ct]] <- df_ct
}
acceleration_oof_df <- do.call(rbind.data.frame, acceleration_oof_df)
rownames(acceleration_oof_df) <- NULL

# Clock performance in controls with out-of-batch predictions
oof_perf <- acceleration_oof_df %>%
        filter(is_ctrl) %>%
        group_by(cell_type) %>%
        summarize(
                r = cor(age, age_pred, use = "complete.obs"),
                r2 = get_r2(age[!is.na(age_pred)], age_pred[!is.na(age_pred)]),
                mae = mean(abs(age - age_pred), na.rm = TRUE),
                n_prots = length(clock_feats[[unique(cell_type)]])
        ) %>%
        arrange(desc(r))
oof_perf

# Agreement with the original accelerations
accel_compare_df <- dplyr::inner_join(
        acceleration_df[, c("sample", "cell_type", "acceleration")],
        acceleration_oof_df[, c("sample", "cell_type", "acceleration", "is_ctrl")],
        by = c("sample", "cell_type"),
        suffix = c("_orig", "_oof")
)
accel_compare_df %>%
        group_by(cell_type, is_ctrl) %>%
        summarize(r = cor(acceleration_orig, acceleration_oof, use = "complete.obs"),
                  .groups = "drop")

p_oof_vs_orig <- ggplot(accel_compare_df,
                        aes(x = acceleration_orig, y = acceleration_oof, color = is_ctrl)) +
        geom_point(size = 0.6, alpha = 0.5) +
        geom_abline(linetype = "dashed") +
        facet_wrap(~cell_type, scales = "free") +
        labs(x = "Acceleration (original)", y = "Acceleration (out-of-batch)",
             color = "Control") +
        theme_classic()
p_oof_vs_orig
ggsave(file.path(meet_dir, "A_oof_vs_original_acceleration.pdf"), p_oof_vs_orig,
       height = 6, width = 8)

# Disease associations with out-of-batch accelerations. Same model as the main
# analysis (adjusted for global acceleration and QAlb), plus age and sex.
acceleration_oof_df$accel_global <- acceleration_oof_df$acceleration[
        match(
                paste(acceleration_oof_df$sample, "Global"),
                paste(acceleration_oof_df$sample, acceleration_oof_df$cell_type)
        )
]
cell_clocks <- setdiff(unique(acceleration_oof_df$cell_type), "Global")

assoc_oof <- assoc_by_clock(
        acceleration_oof_df,
        "diagnosis_detailed + accel_global + QAlb_log2 + age + sex",
        clocks = cell_clocks
) %>%
        filter(grepl("^diagnosis_detailed", term)) %>%
        mutate(disease = gsub("^diagnosis_detailed", "", term)) %>%
        group_by(cell_type) %>%
        mutate(p_adj = p.adjust(p, method = "BH")) %>%
        ungroup()

assoc_orig <- disease_cell_assoc_wholemod_subdiag %>%
        filter(grepl("^diagnosis_detailed", disease)) %>%
        transmute(cell_type,
                  disease = gsub("^diagnosis_detailed", "", disease),
                  estimate_orig = Estimate)

assoc_oof_vs_orig <- dplyr::inner_join(assoc_oof, assoc_orig, by = c("cell_type", "disease"))

p_assoc_oof_vs_orig <- ggplot(assoc_oof_vs_orig,
                              aes(x = estimate_orig, y = estimate, color = cell_type)) +
        geom_abline(linetype = "dashed", color = "grey50") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
        geom_point(size = 2.5) +
        ggrepel::geom_text_repel(aes(label = disease), size = 2.5,
                                 max.overlaps = 20, show.legend = FALSE) +
        labs(x = "Disease effect (years), original",
             y = "Disease effect (years), out-of-batch",
             color = "Cell type") +
        theme_classic()
p_assoc_oof_vs_orig
ggsave(file.path(meet_dir, "A_disease_effects_oof_vs_original.pdf"),
       p_assoc_oof_vs_orig, height = 5, width = 6)

p_hm_oof <- ggplot(assoc_oof %>% filter(disease %in% diag_model),
                   aes(x = disease, y = cell_type, fill = estimate)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_stars(p_adj)), vjust = 0.75, size = 4) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "Acceleration\n(years)") +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_hm_oof
ggsave(file.path(meet_dir, "A_cell_disease_hm_oof.pdf"), p_hm_oof,
       height = 3.5, width = 6)

write.csv(oof_perf, file.path(meet_dir, "A_oof_clock_performance.csv"), row.names = FALSE)
write.csv(acceleration_oof_df, file.path(meet_dir, "A_acceleration_oof.csv"), row.names = FALSE)
write.csv(assoc_oof, file.path(meet_dir, "A_disease_assoc_oof.csv"), row.names = FALSE)

# Wide matrices used in the following sections
accel_oof_wide <- acceleration_oof_df %>%
        select(sample, cell_type, acceleration) %>%
        tidyr::pivot_wider(names_from = cell_type, values_from = acceleration) %>%
        as.data.frame()
rownames(accel_oof_wide) <- accel_oof_wide$sample
accel_oof_wide <- accel_oof_wide[samp_df$sample, -1]

accel_z_wide <- acceleration_oof_df %>%
        select(sample, cell_type, accel_z) %>%
        tidyr::pivot_wider(names_from = cell_type, values_from = accel_z) %>%
        as.data.frame()
rownames(accel_z_wide) <- accel_z_wide$sample
accel_z_wide <- accel_z_wide[samp_df$sample, -1]

################################################################################
# B. Aging-axis alignment                                                      #
################################################################################
# For each protein: age effect in controls (SD/year) and effect of each disease
# vs controls (SD), both adjusted for age/sex. For each cell type, correlate the
# two across the clock proteins: rho > 0 means the disease moves that cell
# type's proteome in the direction of aging ("accelerated aging"), rho < 0 the
# opposite, ~0 an orthogonal (non-aging) disease process. The slope of disease
# effects on age effects reads as "equivalent years of aging".
# Proteins are correlated, so instead of cor.test p-values we compare each
# cell type with random protein sets of the same size (empirical p).

# Standardize every protein using the controls
ctrl_mu <- colMeans(bader_bc_q[samp_df$is_ctrl, ], na.rm = TRUE)
ctrl_sd <- apply(bader_bc_q[samp_df$is_ctrl, ], 2, sd, na.rm = TRUE)
ok_prot <- colSums(!is.na(bader_bc_q[samp_df$is_ctrl, ])) >= 30 & ctrl_sd > 0 & !is.na(ctrl_sd)
P_std <- sweep(sweep(bader_bc_q[, ok_prot], 2, ctrl_mu[ok_prot], "-"), 2, ctrl_sd[ok_prot], "/")

# Coefficient of a term for every column of P, with per-protein NA handling
coef_per_protein <- function(P, design, term){
        out <- rep(NA_real_, ncol(P))
        names(out) <- colnames(P)
        for (j in seq_len(ncol(P))){
                yv <- P[, j]
                ok <- !is.na(yv) & complete.cases(design)
                if (sum(ok) < ncol(design) + 10) next
                fit <- lm.fit(design[ok, , drop = FALSE], yv[ok])
                out[j] <- fit$coefficients[term]
        }
        out
}

print("Computing protein age effects in controls...")
model_matrix_na <- function(f, data){
        model.matrix(f, model.frame(f, data, na.action = na.pass))
}
des_ctrl <- model_matrix_na(~ age + sex, samp_df[samp_df$is_ctrl, ])
beta_age <- coef_per_protein(P_std[samp_df$is_ctrl, ], des_ctrl, "age")

beta_dis <- list()
for (d in diag_model){
        print(sprintf("Computing protein effects of %s...", d))
        idx <- samp_df$is_ctrl | samp_df$diagnosis_detailed == d
        sub <- samp_df[idx, ]
        sub$grp <- as.numeric(sub$diagnosis_detailed == d)
        des_d <- model_matrix_na(~ grp + age + sex, sub)
        beta_dis[[d]] <- coef_per_protein(P_std[idx, ], des_d, "grp")
}
beta_dis <- do.call(cbind, beta_dis)

prot_effects_df <- data.frame(protein = names(beta_age), beta_age = beta_age, beta_dis)
prot_effects_df$cell_type <- prot_info$cell_type[match(prot_effects_df$protein, prot_info$protein)]
write.csv(prot_effects_df, file.path(meet_dir, "B_protein_age_disease_effects.csv"), row.names = FALSE)

alignment_stat <- function(ba, bd){ 
        ok <- !is.na(ba) & !is.na(bd)
        c(rho = cor(ba[ok], bd[ok], method = "spearman"),
          eq_years = unname(coef(lm(bd[ok] ~ 0 + ba[ok]))[1]))
}

n_null <- 5000
set.seed(42)
align_res <- list()
for (ct in cell_clocks){
        prots_ct <- intersect(clock_feats[[ct]], names(beta_age))
        for (d in diag_model){
                obs <- alignment_stat(beta_age[prots_ct], beta_dis[prots_ct, d])
                null_rho <- replicate(n_null, {
                        rp <- sample(names(beta_age), length(prots_ct))
                        alignment_stat(beta_age[rp], beta_dis[rp, d])["rho"]
                })
                align_res[[paste(ct, d)]] <- data.frame(
                        cell_type = ct,
                        disease = d,
                        n_prots = length(prots_ct),
                        rho = obs["rho"],
                        eq_years = obs["eq_years"],
                        null_rho_mean = mean(null_rho, na.rm = TRUE),
                        p_emp = (sum(abs(null_rho - mean(null_rho, na.rm = TRUE)) >=
                                             abs(obs["rho"] - mean(null_rho, na.rm = TRUE)),
                                     na.rm = TRUE) + 1) / (n_null + 1)
                )
        }
}
align_res <- do.call(rbind.data.frame, align_res)
rownames(align_res) <- NULL
align_res$p_adj <- p.adjust(align_res$p_emp, method = "BH")
align_res

p_align <- ggplot(align_res, aes(x = disease, y = cell_type, fill = rho)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_stars(p_adj)), vjust = 0.75, size = 4) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "Spearman rho\n(age vs disease\neffects)") +
        labs(x = NULL, y = NULL,
             subtitle = "Stars: empirical p vs random protein sets (BH)") +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_align
ggsave(file.path(meet_dir, "B_aging_axis_alignment_hm.pdf"), p_align,
       height = 3.5, width = 6.5)

# Scatter for the most aligned pair, as an example
top_pair <- align_res[which.max(abs(align_res$rho)), ]
prots_top <- intersect(clock_feats[[top_pair$cell_type]], names(beta_age))
p_align_ex <- ggplot(
        data.frame(beta_age = beta_age[prots_top] * 10,
                   beta_dis = beta_dis[prots_top, top_pair$disease]),
        aes(x = beta_age, y = beta_dis)) +
        geom_hline(yintercept = 0, color = "grey70") +
        geom_vline(xintercept = 0, color = "grey70") +
        geom_point() +
        geom_smooth(method = "lm") +
        labs(x = "Age effect in controls (SD / decade)",
             y = sprintf("%s vs control (SD)", top_pair$disease),
             title = sprintf("%s clock proteins", top_pair$cell_type)) +
        theme_classic()
p_align_ex
ggsave(file.path(meet_dir, "B_aging_axis_alignment_example.pdf"), p_align_ex,
       height = 4, width = 4.5)
write.csv(align_res, file.path(meet_dir, "B_aging_axis_alignment.csv"), row.names = FALSE)

################################################################################
# C. Pseudo-temporal dynamics from cross-sectional data                        #
################################################################################

# C1. Does the gap between patients and controls widen with age?
# (diagnosis x age interaction; estimates are years of acceleration per decade)
acceleration_oof_df$age_dec <- acceleration_oof_df$age / 10
assoc_age_int <- assoc_by_clock(
        acceleration_oof_df[acceleration_oof_df$diagnosis_detailed %in% c(ctrl_lab, diag_model), ],
        "diagnosis_detailed * age_dec + sex",
        clocks = cell_clocks
) %>%
        filter(grepl(":age_dec$", term)) %>%
        mutate(disease = gsub("^diagnosis_detailed|:age_dec$", "", term),
               p_adj = p.adjust(p, method = "BH"))

p_age_int <- ggplot(assoc_age_int, aes(x = disease, y = cell_type, fill = estimate)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_stars(p_adj)), vjust = 0.75, size = 4) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "Change in\nacceleration\n(years/decade)") +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_age_int
ggsave(file.path(meet_dir, "C1_diagnosis_by_age_interaction_hm.pdf"), p_age_int,
       height = 3.5, width = 6)

p_traj <- ggplot(
        acceleration_oof_df[acceleration_oof_df$diagnosis_detailed %in% c(ctrl_lab, diag_model) &
                                    acceleration_oof_df$cell_type %in% cell_clocks, ],
        aes(x = age, y = acceleration, color = diagnosis_detailed)) +
        geom_hline(yintercept = 0, color = "grey70") +
        geom_smooth(method = "loess", se = FALSE, span = 1) +
        facet_wrap(~cell_type, scales = "free_y") +
        labs(x = "Age", y = "Acceleration (years)", color = NULL) +
        theme_classic()
p_traj
ggsave(file.path(meet_dir, "C1_acceleration_trajectories_by_age.pdf"), p_traj,
       height = 6, width = 9)

# C2. Abnormality-frequency ordering (event-based-model logic): in a disease,
# cell types that are abnormal in more patients are inferred to become
# abnormal earlier. Here: excess fraction of patients with z >= threshold over
# controls. RRMS vs PMS gives a first "stage" contrast. With longitudinal
# Knight ADRC data this could be done with a proper EBM / SuStaIn.
abn_freq <- acceleration_oof_df %>%
        filter(cell_type %in% cell_clocks,
               diagnosis_detailed %in% c(ctrl_lab, diag_model)) %>%
        group_by(cell_type, diagnosis_detailed) %>%
        summarize(frac_abn = mean(accel_z >= z_thr_gap, na.rm = TRUE),
                  n = sum(!is.na(accel_z)), .groups = "drop") %>%
        group_by(cell_type) %>%
        mutate(excess = frac_abn - frac_abn[diagnosis_detailed == ctrl_lab]) %>%
        ungroup() %>%
        filter(diagnosis_detailed != ctrl_lab)

p_abn <- ggplot(abn_freq,
                aes(x = excess,
                    y = reorder(paste0(cell_type, "___", diagnosis_detailed), excess))) +
        geom_vline(xintercept = 0, color = "grey70") +
        geom_col(fill = "#3B6E8C") +
        facet_wrap(~diagnosis_detailed, scales = "free_y") +
        scale_y_discrete(labels = function(x) sub("___.*", "", x)) +
        labs(x = sprintf("Excess fraction of patients with z >= %s vs controls", z_thr_gap),
             y = NULL) +
        theme_classic()
p_abn
ggsave(file.path(meet_dir, "C2_abnormality_frequency_ordering.pdf"), p_abn,
       height = 6, width = 9)
write.csv(abn_freq, file.path(meet_dir, "C2_abnormality_frequency.csv"), row.names = FALSE)

# C3. Clinical / progression variables available in the metadata. We scan every
# numeric column (excluding identifiers and the variables already modelled)
# and test its association with each cell-type acceleration: across all
# samples (adjusted for age, sex, diagnosis) and within each diagnosis.
# Check `clin_vars` and drop/add columns as needed (e.g. disease duration,
# EDSS, MMSE, CSF cell count, lactate, IgG index, NfL...).
excl_cols <- c("Age", "Sample_preparation_batch", "Sample_plate", "ID",
               "Run_timsTOFMeasurement", "sample")
clin_vars <- colnames(bader_mdat)[sapply(bader_mdat, function(v){
        is.numeric(v) && sum(!is.na(v)) >= 50 && length(unique(na.omit(v))) > 5
})]
clin_vars <- setdiff(clin_vars, excl_cols)
clin_vars <- clin_vars[!grepl("batch|plate|run|order|position|index_in|^X$|id$",
                              clin_vars, ignore.case = TRUE)]
clin_vars

# Progression-like variables, flagged for the meeting
prog_vars <- clin_vars[grepl("durat|onset|edss|mmse|moca|cdr|stage|progr|sever|follow|nfl",
                             clin_vars, ignore.case = TRUE)]
prog_vars

clin_res <- list()
for (v in clin_vars){
        for (ct in cell_clocks){
                df_v <- acceleration_oof_df[acceleration_oof_df$cell_type == ct, ]
                df_v$clin <- as.numeric(scale(bader_mdat[match(df_v$sample, bader_mdat$Run_timsTOFMeasurement), v]))
                # All samples
                fit <- try(lm(scale(acceleration) ~ clin + age + sex + diagnosis_detailed, data = df_v),
                           silent = TRUE)
                if (!inherits(fit, "try-error") && "clin" %in% rownames(coef(summary(fit)))){
                        cf <- coef(summary(fit))["clin", ]
                        clin_res[[paste(v, ct, "all")]] <- data.frame(
                                variable = v, cell_type = ct, group = "All (adj. diagnosis)",
                                n = nobs(fit), estimate = cf[1], p = cf[4])
                }
                # Within each diagnosis
                for (d in c(ctrl_lab, diag_model)){
                        df_vd <- df_v[df_v$diagnosis_detailed == d, ]
                        if (sum(!is.na(df_vd$clin)) < 20) next
                        fit <- try(lm(scale(acceleration) ~ clin + age + sex, data = df_vd),
                                   silent = TRUE)
                        if (inherits(fit, "try-error") || !"clin" %in% rownames(coef(summary(fit)))) next
                        cf <- coef(summary(fit))["clin", ]
                        clin_res[[paste(v, ct, d)]] <- data.frame(
                                variable = v, cell_type = ct, group = d,
                                n = nobs(fit), estimate = cf[1], p = cf[4])
                }
        }
}
if (length(clin_res) > 0){
        clin_res <- do.call(rbind.data.frame, clin_res)
        rownames(clin_res) <- NULL
        clin_res$p_adj <- p.adjust(clin_res$p, method = "BH")
        clin_res$is_progression_var <- clin_res$variable %in% prog_vars
        print(head(clin_res[order(clin_res$p), ], 30))
        write.csv(clin_res, file.path(meet_dir, "C3_clinical_variable_associations.csv"),
                  row.names = FALSE)
        
        p_clin <- ggplot(clin_res[clin_res$group == "All (adj. diagnosis)", ],
                         aes(x = variable, y = cell_type, fill = estimate)) +
                geom_tile(color = "white", linewidth = 0.6) +
                geom_text(aes(label = sig_stars(p_adj)), vjust = 0.75, size = 4) +
                scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                     midpoint = 0, name = "Std. beta") +
                labs(x = NULL, y = NULL) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        p_clin
        ggsave(file.path(meet_dir, "C3_clinical_variable_associations_hm.pdf"), p_clin,
               height = 4, width = max(5, 0.4 * length(clin_vars) + 3))
}

################################################################################
# D. Age-gappers                                                               #
################################################################################
# Extreme agers per cell type (z >= threshold, relative to the controls'
# acceleration SD). D1: are patients enriched among extreme agers of specific
# cell types (OR adjusted for age and sex)? D2: dose-response with the number
# of cell types that are aged in the same person ("multi-cell agers"),
# analogous to multi-organ aging in plasma.

gap_df <- samp_df
gap_df <- cbind(gap_df, (accel_z_wide[, cell_clocks] >= z_thr_gap) * 1)
gap_df$n_aged <- rowSums(accel_z_wide[, cell_clocks] >= z_thr_gap, na.rm = TRUE)
gap_df$n_youthful <- rowSums(accel_z_wide[, cell_clocks] <= -z_thr_gap, na.rm = TRUE)

# D1. Enrichment of each diagnosis among extreme agers of each cell type
gap_or <- list()
for (d in diag_model){
        sub <- gap_df[gap_df$is_ctrl | gap_df$diagnosis_detailed == d, ]
        sub$case <- as.numeric(sub$diagnosis_detailed == d)
        for (ct in cell_clocks){
                sub$ager <- sub[, ct]
                if (length(unique(na.omit(sub$ager))) < 2) next
                fit <- glm(case ~ ager + age + sex, data = sub, family = binomial)
                cf <- coef(summary(fit))
                if (!"ager" %in% rownames(cf)) next
                gap_or[[paste(d, ct)]] <- data.frame(
                        disease = d, cell_type = ct,
                        n_agers_case = sum(sub$ager[sub$case == 1], na.rm = TRUE),
                        n_agers_ctrl = sum(sub$ager[sub$case == 0], na.rm = TRUE),
                        log_or = cf["ager", 1], se = cf["ager", 2], p = cf["ager", 4]
                )
        }
}
gap_or <- do.call(rbind.data.frame, gap_or)
rownames(gap_or) <- NULL
gap_or$p_adj <- p.adjust(gap_or$p, method = "BH")
gap_or$OR <- exp(gap_or$log_or)
gap_or

p_gap_or <- ggplot(gap_or, aes(x = disease, y = cell_type, fill = log2(OR))) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = sig_stars(p_adj)), vjust = 0.75, size = 4) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "log2 OR\n(extreme ager)",
                             limits = c(-4, 4), oob = scales::squish) +
        labs(x = NULL, y = NULL,
             subtitle = sprintf("Extreme ager: acceleration z >= %s", z_thr_gap)) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_gap_or
ggsave(file.path(meet_dir, "D1_extreme_ager_enrichment_hm.pdf"), p_gap_or,
       height = 3.5, width = 6)
write.csv(gap_or, file.path(meet_dir, "D1_extreme_ager_enrichment.csv"), row.names = FALSE)

# D2. Multi-cell agers: prevalence of each diagnosis vs number of aged cell types
gap_df$n_aged_cat <- factor(pmin(gap_df$n_aged, 3), levels = 0:3,
                            labels = c("0", "1", "2", "3+"))
table(gap_df$n_aged_cat, gap_df$diagnosis_detailed)

multi_res <- list()
for (d in diag_model){
        sub <- gap_df[gap_df$is_ctrl | gap_df$diagnosis_detailed == d, ]
        sub$case <- as.numeric(sub$diagnosis_detailed == d)
        fit <- glm(case ~ n_aged + n_youthful + age + sex, data = sub, family = binomial)
        cf <- coef(summary(fit))
        multi_res[[d]] <- data.frame(
                disease = d,
                term = c("n_aged", "n_youthful"),
                OR_per_cell_type = exp(cf[c("n_aged", "n_youthful"), 1]),
                p = cf[c("n_aged", "n_youthful"), 4]
        )
}
multi_res <- do.call(rbind.data.frame, multi_res)
rownames(multi_res) <- NULL
multi_res$p_adj <- p.adjust(multi_res$p, method = "BH")
multi_res

multi_prev <- gap_df %>%
        filter(diagnosis_detailed %in% c(ctrl_lab, diag_model)) %>%
        group_by(n_aged_cat) %>%
        mutate(n_cat = n()) %>%
        group_by(n_aged_cat, diagnosis_detailed, n_cat) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(prop = n / n_cat)

p_multi <- ggplot(multi_prev, aes(x = n_aged_cat, y = prop, fill = diagnosis_detailed)) +
        geom_col() +
        geom_text(data = distinct(multi_prev, n_aged_cat, n_cat),
                  aes(x = n_aged_cat, y = 1.03, label = paste0("n=", n_cat)),
                  inherit.aes = FALSE, size = 3) +
        labs(x = "Number of aged cell types", y = "Proportion", fill = NULL) +
        theme_classic()
p_multi
ggsave(file.path(meet_dir, "D2_multi_cell_agers_composition.pdf"), p_multi,
       height = 4, width = 5.5)
write.csv(multi_res, file.path(meet_dir, "D2_multi_cell_agers_OR.csv"), row.names = FALSE)

# Co-occurrence of extreme aging across cell types (in controls only, to avoid
# disease-driven correlation): are there "ageotypes"?
co_age <- cor(accel_oof_wide[samp_df$is_ctrl, cell_clocks], use = "pairwise.complete.obs")
co_age
pdf(file.path(meet_dir, "D3_cell_type_acceleration_correlation_controls.pdf"),
    height = 5, width = 5.5)
heatmap(co_age, symm = TRUE, scale = "none",
        col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(51),
        zlim = c(-1, 1), margins = c(8, 8))
dev.off()

################################################################################
# E. Within age-gappers, do clock proteins separate disease from healthy      #
#    agers?                                                                    #
################################################################################
# Cross-sectional analogue of the UKBB analysis. Among individuals with an
# accelerated clock for cell type X (z >= soft threshold), classify patients vs
# controls using:
#   1. Acceleration of that clock + age + sex   (the age gap alone)
#   2. The clock proteins (+ age, sex unpenalized)
#   3. Random protein sets of the same size     (does the cell-type signature
#      matter, or is any set of CSF proteins enough?)
# Outcomes: any disease vs control, and each diagnosis vs control.

n_rand_sets <- 10
min_class_n <- 10

gapper_clf <- list()
set.seed(2026)
for (ct in cell_clocks){
        print(sprintf("Age-gapper classification: %s...", ct))
        is_gapper <- !is.na(accel_z_wide[, ct]) & accel_z_wide[, ct] >= z_thr_gap_soft
        prots_ct <- intersect(clock_feats[[ct]], colnames(bader_bc))
        for (outc in c("Any disease", diag_model)){
                if (outc == "Any disease"){
                        idx <- is_gapper
                        yv <- as.numeric(!samp_df$is_ctrl[idx])
                }else{
                        idx <- is_gapper & (samp_df$is_ctrl | samp_df$diagnosis_detailed == outc)
                        yv <- as.numeric(samp_df$diagnosis_detailed[idx] == outc)
                }
                if (min(table(factor(yv, levels = 0:1))) < min_class_n) next
                cov_df <- data.frame(age = samp_df$age[idx], sex = samp_df$sex[idx])
                keep <- complete.cases(cov_df)
                idx_w <- which(idx)[keep]
                yv <- yv[keep]
                cov_df <- droplevels(cov_df[keep, , drop = FALSE])
                
                res_gap <- cv_auc(yv, covars = cbind(cov_df, accel = accel_oof_wide[idx_w, ct]))
                res_prot <- cv_auc(yv, X = bader_bc_q[idx_w, prots_ct, drop = FALSE], covars = cov_df)
                res_rand <- sapply(seq_len(n_rand_sets), function(i){
                        rp <- sample(setdiff(colnames(P_std), prots_ct), length(prots_ct))
                        cv_auc(yv, X = bader_bc_q[idx_w, rp, drop = FALSE], covars = cov_df, n_rep = 2)["auc"]
                })
                gapper_clf[[paste(ct, outc)]] <- data.frame(
                        cell_type = ct, outcome = outc,
                        n_case = sum(yv == 1), n_ctrl = sum(yv == 0),
                        auc_age_gap = res_gap["auc"],
                        auc_clock_prots = res_prot["auc"],
                        auc_clock_prots_sd = res_prot["auc_sd"],
                        auc_random_prots = mean(res_rand),
                        auc_random_prots_sd = sd(res_rand)
                )
        }
}
gapper_clf <- do.call(rbind.data.frame, gapper_clf)
rownames(gapper_clf) <- NULL
gapper_clf

gapper_clf_long <- gapper_clf %>%
        select(cell_type, outcome, auc_age_gap, auc_clock_prots, auc_random_prots) %>%
        tidyr::pivot_longer(cols = starts_with("auc_"), names_to = "model", values_to = "auc") %>%
        mutate(model = factor(gsub("auc_", "", model),
                              levels = c("age_gap", "random_prots", "clock_prots")))

p_gapper <- ggplot(gapper_clf_long, aes(x = cell_type, y = auc, fill = model)) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
        geom_col(position = position_dodge(0.8), width = 0.75) +
        coord_flip(ylim = c(0.4, 1)) +
        facet_wrap(~outcome) +
        scale_fill_manual(values = c(age_gap = "grey70", random_prots = "#9DB4C0",
                                     clock_prots = "#3B6E8C")) +
        labs(x = NULL, y = "Cross-validated AUC", fill = NULL,
             subtitle = sprintf("Individuals with acceleration z >= %s in each cell type",
                                z_thr_gap_soft)) +
        theme_classic()
p_gapper
ggsave(file.path(meet_dir, "E_age_gapper_disease_classification.pdf"), p_gapper,
       height = 6, width = 10)
write.csv(gapper_clf, file.path(meet_dir, "E_age_gapper_disease_classification.csv"),
          row.names = FALSE)

# Proteins driving the separation in the best case (fit on all gappers)
best_case <- gapper_clf[which.max(gapper_clf$auc_clock_prots - gapper_clf$auc_random_prots), ]
best_case
is_gapper <- !is.na(accel_z_wide[, best_case$cell_type]) &
        accel_z_wide[, best_case$cell_type] >= z_thr_gap_soft
if (best_case$outcome == "Any disease"){
        idx <- is_gapper
        yv <- as.numeric(!samp_df$is_ctrl[idx])
}else{
        idx <- is_gapper & (samp_df$is_ctrl | samp_df$diagnosis_detailed == best_case$outcome)
        yv <- as.numeric(samp_df$diagnosis_detailed[idx] == best_case$outcome)
}
prots_best <- intersect(clock_feats[[best_case$cell_type]], colnames(bader_bc))
X_best <- impute_train_test(bader_bc_q[idx, prots_best, drop = FALSE],
                            bader_bc_q[idx, prots_best, drop = FALSE])$train
fit_best <- cv.glmnet(scale(as.matrix(X_best)), yv, family = "binomial", alpha = 0.5)
coef_best <- as.matrix(coef(fit_best, s = "lambda.min"))
coef_best <- data.frame(protein = rownames(coef_best), coef = coef_best[, 1])
coef_best <- coef_best[coef_best$protein != "(Intercept)" & coef_best$coef != 0, ]
coef_best$beta_age <- beta_age[coef_best$protein]
coef_best <- coef_best[order(-abs(coef_best$coef)), ]
head(coef_best, 20)
write.csv(coef_best, file.path(meet_dir, "E_best_case_signature.csv"), row.names = FALSE)

################################################################################
# F. Added value of the cell-type acceleration profile for diagnosis           #
################################################################################
# Each diagnosis vs controls, CV AUC of logistic models with:
#   age + sex  |  + global acceleration  |  + all cell-type accelerations
diag_clf <- list()
set.seed(7)
for (d in diag_model){
        idx <- samp_df$is_ctrl | samp_df$diagnosis_detailed == d
        base <- data.frame(age = samp_df$age[idx], sex = samp_df$sex[idx])
        A <- accel_oof_wide[idx, , drop = FALSE]
        colnames(A) <- make.names(colnames(A))
        keep <- complete.cases(base) & complete.cases(A)
        yv <- as.numeric(samp_df$diagnosis_detailed[idx] == d)[keep]
        base <- droplevels(base[keep, , drop = FALSE])
        A <- A[keep, , drop = FALSE]
        diag_clf[[d]] <- data.frame(
                disease = d,
                model = c("Age + sex", "+ global clock", "+ cell-type clocks",
                          "+ global + cell-type clocks"),
                auc = c(
                        cv_auc(yv, covars = base)["auc"],
                        cv_auc(yv, covars = cbind(base, Global = A$Global))["auc"],
                        cv_auc(yv, covars = cbind(base, A[, make.names(cell_clocks)]))["auc"],
                        cv_auc(yv, covars = cbind(base, A))["auc"]
                )
        )
}
diag_clf <- do.call(rbind.data.frame, diag_clf)
rownames(diag_clf) <- NULL
diag_clf$model <- factor(diag_clf$model, levels = unique(diag_clf$model))
diag_clf

p_diag_clf <- ggplot(diag_clf, aes(x = disease, y = auc, fill = model)) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
        geom_col(position = position_dodge(0.8), width = 0.75) +
        coord_flip(ylim = c(0.4, 1)) +
        scale_fill_manual(values = c("grey75", "#9DB4C0", "#3B6E8C", "#1B3A4B")) +
        labs(x = NULL, y = "Cross-validated AUC (vs controls)", fill = NULL) +
        theme_classic()
p_diag_clf
ggsave(file.path(meet_dir, "F_diagnosis_classification_from_clocks.pdf"), p_diag_clf,
       height = 4, width = 6)
write.csv(diag_clf, file.path(meet_dir, "F_diagnosis_classification_from_clocks.csv"),
          row.names = FALSE)

################################################################################
# G. Sex differences                                                           #
################################################################################
assoc_sex_int <- assoc_by_clock(
        acceleration_oof_df[acceleration_oof_df$diagnosis_detailed %in% c(ctrl_lab, diag_model), ],
        "diagnosis_detailed * sex + age",
        clocks = c(cell_clocks, "Global")
)
assoc_sex_main <- assoc_sex_int %>%
        filter(grepl("^sex", term)) %>%
        mutate(p_adj = p.adjust(p, method = "BH"))
assoc_sex_int <- assoc_sex_int %>%
        filter(grepl(":sex", term)) %>%
        mutate(disease = gsub("^diagnosis_detailed|:sex.*$", "", term),
               p_adj = p.adjust(p, method = "BH"))
assoc_sex_main
assoc_sex_int[order(assoc_sex_int$p), ]
write.csv(rbind(assoc_sex_main[, colnames(assoc_sex_main)],
                assoc_sex_int[, colnames(assoc_sex_main)]),
          file.path(meet_dir, "G_sex_effects.csv"), row.names = FALSE)

################################################################################
# H. Transferability to SomaScan 7k (Knight ADRC)                              #
################################################################################
# How many proteins of each cell-type clock are measured on the SomaScan 7k
# panel, and how well do clocks perform when restricted to Bader proteins that
# are also on SomaScan? (Same out-of-batch framework, fixed alpha = 0.5 for both
# the full and the restricted candidate sets, so they are comparable.)
soma7k <- as.data.frame(readxl::read_xlsx(somascan7k_f))
soma_uniprot_col <- grep("uniprot", colnames(soma7k), ignore.case = TRUE, value = TRUE)[1]
soma_uniprot <- unique(unlist(strsplit(as.character(soma7k[, soma_uniprot_col]), "[;|, ]+")))
soma_uniprot <- gsub("-.*", "", soma_uniprot[!is.na(soma_uniprot) & soma_uniprot != ""])

bader_on_soma <- sapply(
        strsplit(prot_info$protein, ";"),
        function(v) any(gsub("-.*", "", v) %in% soma_uniprot)
)
names(bader_on_soma) <- prot_info$protein
mean(bader_on_soma)

soma_cov <- list()
set.seed(99)
for (ct in cell_clocks){
        print(sprintf("SomaScan-restricted clock: %s...", ct))
        cand <- intersect(prot_info$protein[prot_info$cell_type %in% ct], colnames(bader_bc))
        cand_soma <- cand[bader_on_soma[cand]]
        clock <- intersect(clock_feats[[ct]], names(bader_on_soma))
        perf_set <- function(prots){
                if (length(prots) < 3) return(NA_real_)
                pr <- lobo_clock_predict(prots, 0.5, bader_bc, samp_df$age,
                                         samp_df$batch, is_train_ctrl)
                cor(pr[is_train_ctrl], samp_df$age[is_train_ctrl], use = "complete.obs")
        }
        soma_cov[[ct]] <- data.frame(
                cell_type = ct,
                n_candidates = length(cand),
                n_candidates_soma = length(cand_soma),
                n_clock = length(clock),
                n_clock_soma = sum(bader_on_soma[clock]),
                r_all_candidates = perf_set(cand),
                r_soma_candidates = perf_set(cand_soma)
        )
}
soma_cov <- do.call(rbind.data.frame, soma_cov)
rownames(soma_cov) <- NULL
soma_cov$frac_clock_soma <- soma_cov$n_clock_soma / soma_cov$n_clock
soma_cov

p_soma <- soma_cov %>%
        mutate(label = sprintf("%s (%d/%d)", cell_type, n_candidates_soma, n_candidates)) %>%
        select(label, r_all_candidates, r_soma_candidates) %>%
        tidyr::pivot_longer(-label, names_to = "set", values_to = "r") %>%
        mutate(set = ifelse(set == "r_all_candidates", "All Bader markers",
                            "Markers also on SomaScan 7k")) %>%
        ggplot(aes(x = reorder(label, r), y = r, fill = set)) +
        geom_col(position = position_dodge(0.8), width = 0.75) +
        coord_flip() +
        scale_fill_manual(values = c("#9DB4C0", "#3B6E8C")) +
        labs(x = NULL, y = "Out-of-batch r (controls)", fill = NULL,
             subtitle = "In brackets: markers on SomaScan / markers in Bader") +
        theme_classic()
p_soma
ggsave(file.path(meet_dir, "H_somascan_transferability.pdf"), p_soma, height = 4, width = 6.5)
write.csv(soma_cov, file.path(meet_dir, "H_somascan_transferability.csv"), row.names = FALSE)

################################################################################
# I. Longitudinal samples in the Bader dataset (if any)                        #
################################################################################
# The clocks were trained on the cross-sectional discovery set. If the full
# supplementary table includes repeated samples of the same individuals, we can
# get a first look at within-person change of cell-type acceleration (the core
# of the "dynamics" paper, to be done properly in Knight ADRC).
table(mdat_all$Experiment)
rep_ids <- names(which(table(mdat_all$ID) > 1))
length(rep_ids)

if (length(rep_ids) > 0){
        long_mdat <- mdat_all[mdat_all$ID %in% rep_ids & !grepl("Pool", mdat_all$Sex), ]
        long_mdat$sample_mn <- make.names(long_mdat$Run_timsTOFMeasurement)
        table(long_mdat$Experiment, long_mdat$Diagnosis_group)
        
        bader_full <- plotUtils::read_table_fast(bader_full_f, row.names = 1)
        rownames(bader_full) <- bader_full$protein
        bader_full <- bader_full[, colnames(bader_full) %in% long_mdat$sample_mn]
        long_mdat <- long_mdat[long_mdat$sample_mn %in% colnames(bader_full), ]
        bader_full <- log2(t(as.matrix(bader_full[, long_mdat$sample_mn])))
        mode(bader_full) <- "numeric"
        
        # Batch-normalize with the cross-sectional controls (unseen batches are
        # left as they are, see batch_normalize)
        common_prots <- intersect(colnames(bader_full), colnames(bader_complete))
        long_bc <- batch_normalize_keep_order(
                bader_complete[samp_df$is_ctrl, common_prots],
                samp_df$batch[samp_df$is_ctrl],
                bader_full[, common_prots],
                long_mdat$Sample_preparation_batch
        )$test
        
        long_accel <- list()
        for (ct in cell_clocks){
                prots <- intersect(clock_feats[[ct]], common_prots)
                tr <- is_train_ctrl & !is.na(samp_df$age)
                res <- qalb_residualize(bader_bc[tr, prots], long_bc[, prots])
                imp <- impute_train_test(res$train, res$test)
                fit <- cv.glmnet(as.matrix(imp$train), samp_df$age[tr],
                                 alpha = clock_alphas[[ct]], nfolds = 5)
                pr <- predict(fit, newx = as.matrix(imp$test), s = "lambda.min")[, 1]
                # Age-bias correction from the out-of-batch control predictions
                ctrl_ct <- acceleration_oof_df[acceleration_oof_df$cell_type == ct &
                                                       acceleration_oof_df$is_ctrl, ]
                bias_fit <- lm(raw_accel ~ age, data = ctrl_ct)
                long_accel[[ct]] <- data.frame(
                        ID = long_mdat$ID,
                        sample = long_mdat$Run_timsTOFMeasurement,
                        diagnosis = long_mdat$Diagnosis_group,
                        age = long_mdat$Age,
                        cell_type = ct,
                        acceleration = (pr - long_mdat$Age) -
                                predict(bias_fit, data.frame(age = long_mdat$Age))
                )
        }
        long_accel <- do.call(rbind.data.frame, long_accel)
        long_accel <- long_accel %>%
                group_by(ID) %>%
                mutate(time = age - min(age, na.rm = TRUE),
                       age_baseline = min(age, na.rm = TRUE)) %>%
                ungroup()
        long_accel$diagnosis <- relevel(factor(long_accel$diagnosis), ref = ctrl_lab)
        
        # Within-person change of acceleration per year, by diagnosis
        long_res <- list()
        for (ct in cell_clocks){
                fit <- try(lmerTest::lmer(
                        acceleration ~ time * diagnosis + age_baseline + (1 | ID),
                        data = long_accel[long_accel$cell_type == ct, ]
                ), silent = TRUE)
                if (inherits(fit, "try-error")) next
                cf <- as.data.frame(coef(summary(fit)))
                cf$term <- rownames(cf)
                cf$cell_type <- ct
                long_res[[ct]] <- cf[grepl("time", cf$term), ]
        }
        long_res <- do.call(rbind.data.frame, long_res)
        long_res
        
        p_long <- ggplot(long_accel, aes(x = age, y = acceleration, group = ID,
                                         color = diagnosis)) +
                geom_line(alpha = 0.5) +
                geom_point(size = 0.8) +
                facet_wrap(~cell_type, scales = "free_y") +
                theme_classic()
        p_long
        ggsave(file.path(meet_dir, "I_longitudinal_acceleration.pdf"), p_long,
               height = 6, width = 9)
        write.csv(long_res, file.path(meet_dir, "I_longitudinal_acceleration_lmm.csv"))
}

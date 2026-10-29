library(caret)


data_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/UKBB"

DGE_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/scRNAseq/DGE"

prot_f <- "proteomics_instance_0.tsv"
mdat_f <- "metadata_filt.csv"

eid_healthy_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/UKBB/healthy_eids.csv"

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

prot <- read.csv(file.path(data_dir, prot_f), sep = '\t')
mdat <- read.csv(file.path(data_dir, mdat_f))

healthy_eids <- read.csv(eid_healthy_f)

colnames(prot)

rownames(prot) <- prot$eid
prot <- prot[, 2:ncol(prot)]

#library(biomaRt)

#mart <- useEnsembl(
#        biomart = "genes",
#        dataset = "hsapiens_gene_ensembl"
#)

#mapping <- getBM(
#        attributes = c("hgnc_symbol", "ensembl_gene_id"),
#        filters = "hgnc_symbol",
#        values = colnames(prot),
#        mart = mart
#)

#write.csv(mapping, file = file.path(data_dir, "ukbb_prot_mapping.csv"))

mapping <- read.csv(file.path(data_dir, "ukbb_prot_mapping.csv"))

na_prop_prot <- apply(prot, 2, function(x) sum(is.na(x))/length(x))
sum(na_prop_prot <= .2)
sum(na_prop_prot >= .2)
na_prop_samp <- apply(prot, 1, function(x) sum(is.na(x))/length(x))
sum(na_prop_samp <= .2)
sum(na_prop_samp >= .8)

# Remove stuff with lots of NAs
prot <- prot[na_prop_samp <= .2, na_prop_prot <= .2]

prot_info <- data.frame(protein = colnames(prot),
                        hgnc_symbol = mapping$hgnc_symbol[match(toupper(colnames(prot)), mapping$hgnc_symbol)],
                        ensembl_gene_id = mapping$ensembl_gene_id[match(toupper(colnames(prot)), mapping$hgnc_symbol)])

# Keep only healthy individuals
prot <- prot[rownames(prot) %in% as.character(healthy_eids$eid), ]

mdat <- mdat[match(rownames(prot), mdat$f.eid), ]

colnames(mdat) <- gsub("f.21003.0.0", "Age", colnames(mdat))


library(rsample)

y <- mdat$Age
#batch <- mdat$Sample_preparation_batch

X <- prot


stopifnot(nrow(X) == length(y))


set.seed(1)
folds <- vfold_cv(mdat, v = 10)


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

get_r2 <- function(y, yhat){
        r2 <- 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
        return(r2)
}

#batches <- unique(batch)

library(glmnet)




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
        selected_alphas <- numeric(nrow(folds))
        for (i in seq_len(nrow(folds))) {
                print(sprintf("Fold %s", i))
                
                split_i <- folds$splits[[i]]
                
                train_df <- analysis(split_i)
                test_df  <- assessment(split_i)
                
                # ----------------------------
                # extract matrices
                # ----------------------------
                X_train <- as.matrix(X_cell[as.character(train_df$f.eid), ])
                X_test  <- as.matrix(X_cell[as.character(test_df$f.eid), ])
                
                y_train <- train_df$Age
                y_test  <- test_df$Age
                
                #batch_train <- train_df$Sample_preparation_batch
                #batch_test  <- test_df$Sample_preparation_batch
                
                # ----------------------------
                # batch correction (train only fit)
                # ----------------------------
                #bc <- batch_normalize(X_train, batch_train, X_test, batch_test)
                
                # ----------------------------
                # imputation
                # ----------------------------
                imp <- impute_train_test(X_train, X_test, method = "mean")
                
                Xtr <- imp$train
                Xte <- imp$test
                
                ##########################################################################
                # GLMNET
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
                
                r2_glm_cv <- c(r2_glm_cv, get_r2(y_test, pred_glm))
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
                
                #dtrain <- xgb.DMatrix(as.matrix(Xtr), label = y_train)
                #dtest  <- xgb.DMatrix(as.matrix(Xte), label = y_test)
                
                #xgb_fit <- xgb.train(params, dtrain, nrounds = 200, verbose = 0)
                
                #pred_xgb <- predict(xgb_fit, dtest)
                
                #r2_xgb_cv <- c(r2_xgb_cv, get_r2(y_test, pred_xgb))
                #results_xgb <- c(results_xgb, cor(pred_xgb, y_test, use="complete.obs"))
        }
        
        final_alpha <- as.numeric(
                names(
                        which.max(
                                table(selected_alphas)
                        )
                )
        )
        
        
        
        mean(results_glm)
        #mean(results_xgb)
        mean(r2_glm_cv)
        #mean(r2_xgb_cv)
        
        prot_stats_glm_cv$times_selected <- prot_stats_glm_cv$times_selected / length(seq_len(nrow(folds)))
        
        sel_prots <- prot_stats_glm_cv$protein[prot_stats_glm_cv$times_selected >= .5]
        
        #X_cell_bc <- batch_normalize(
        #        X_cell,
        #        batch,
        #        X_cell,
        #        batch)
        #X_cell_bc <- X_cell_bc$train
        
        # 2. imputation (train only)
        X_cell_imp <- impute_train_test(X_cell, X_cell, method = "mean")
        
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
                fold =  seq_len(nrow(folds)),
                glm_cor = results_glm,
                glm_r2 = r2_glm_cv
        )
        cv_res_list[[ct]] <- cv_df
        final_mod_perf <- data.frame(cell_type = ct, glm_cor = glm_cor_final, glm_r2 = glm_r2_final)
        final_mod_res[[ct]] <- final_mod_perf
        glm_list[[ct]] <- final_cvfit
}

final_mod_res <- do.call(rbind.data.frame, final_mod_res)

cv_res_list <- do.call(rbind.data.frame, cv_res_list)

library(dplyr)
cv_res_mean <- cv_res_list %>%
        group_by(cell_type) %>%
        dplyr::summarize(mean_glm_cor = mean(glm_cor),
                         mean_glm_r2 = mean(glm_r2))

final_mod_res
cv_res_mean

if (!require("biomaRt", quietly = TRUE)) {
        BiocManager::install("biomaRt", update = FALSE)
}
if (!require("readxl", quietly = TRUE)) {
        install.packages("readxl")
}
library(biomaRt)
library(dplyr)
library(tidyr)
library(tibble)


map_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_diversity/emory_div_prot_mapping.csv"
bader_timstof_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/DA-F08.4_-SEC-pass_v06Sc_ion_LibPGQVal1perc_precdLFQv03s01_prot_preprSc04.tsv"
bader_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/1-s2.0-S0092867426001042-mmc2.xlsx"

#emory_div_prot <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_diversity/TMT_Processed_Data/TMT_NormalizedAbundances_NoCorrection.csv", row.names = 1)
emory_div_prot <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_diversity/TMT_Processed_Data/cleanDat.postTAMPORBatchCorrection.csv", row.names = 1)

plot(density(as.matrix(emory_div_prot), na.rm = TRUE))

emory_div_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_diversity/TMT_Processed_Data/DiversityCSF_Traits.csv", row.names = 1)

out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/merged_csf_prots"

plotUtils::create_dir_if_not(out_dir)

mapping <- read.csv(map_file, row.names = 1)
bader_timstof <- plotUtils::read_table_fast(bader_timstof_f, row.names = 1)
bader_mdat <- data.frame(readxl::read_xlsx(bader_mdat_f, sheet = 4))



emory_div_mart_f <- file.path(out_dir, "emory_div_mapping.csv")
if (!file.exists(emory_div_mart_f)){
        mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        
        map <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = mapping$Accession,
                mart = mart
        )
        write.csv(map, file = emory_div_mart_f)
}else{
        map <- read.csv(emory_div_mart_f)
}


# Identify which are the columns that are character
char_cols <- c()
for (i in 1:ncol(bader_timstof)){
        if (class(bader_timstof[, i]) == "character"){
                char_cols <- c(char_cols, i)
        }
}

bader_prot_info <- bader_timstof[, char_cols]
length(intersect(mapping$Accession, bader_prot_info$protein))
length(mapping$Accession)
length(bader_prot_info$protein)


missing_prots <- mapping$Accession[!mapping$Accession %in% bader_prot_info$protein]


emory_div_missing_mart_f <- file.path(out_dir, "emory_div_missing_mapping.csv")
if (!file.exists(emory_div_missing_mart_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        map_missing <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = missing_prots,
                mart = mart
        )
        write.csv(map_missing, file = emory_div_missing_mart_f)
}else{
        map_missing <- read.csv(emory_div_missing_mart_f)
}

# Half of the proteins that don't map don't do so because in the Bader
# dataset there are several proteins that have suffixes such as -2. Remove
# them for these
sum(map_missing$uniprotswissprot %in% gsub("\\-.*", "", bader_prot_info$protein))
any(map_missing$uniprotswissprot %in% gsub("\\-.*", "", bader_prot_info$protein))
any(map_missing$uniprotsptrembl %in% bader_prot_info$protein)


recovered_prots <- map_missing$uniprotswissprot[map_missing$uniprotswissprot %in% gsub("\\-.*", "", bader_prot_info$protein)]


recovered_prots_isoforms <- sort(
        bader_prot_info$protein[gsub("\\-.*", "", bader_prot_info$protein) %in% recovered_prots]
)
length(unique(gsub("\\-.*", "", recovered_prots_isoforms)))
length(gsub("\\-.*", "", recovered_prots_isoforms))


# Keep only sample columns of Bader dataset and filter metadata
bader_timstof <- bader_timstof[, !1:ncol(bader_timstof) %in% char_cols]
rownames(bader_timstof) <- bader_prot_info$protein

bader_timstof <- bader_timstof[, colnames(bader_timstof) %in% make.names(bader_mdat$Run_timsTOFMeasurement)]

bader_mdat <- bader_mdat[match(colnames(bader_timstof), make.names(make.names(bader_mdat$Run_timsTOFMeasurement))), ]
all(colnames(bader_timstof) == make.names(bader_mdat$Run_timsTOFMeasurement))
keep_vec <- bader_mdat$Sample_plate %in% paste0("plate", 1:56) &
        bader_mdat$Experiment == "cross-sectional_discovery" &
        !grepl("Pool", bader_mdat$Sex)
bader_timstof <- bader_timstof[, keep_vec]
bader_mdat <- bader_mdat[keep_vec, ]


plot(density(as.matrix(bader_timstof), na.rm = TRUE))


bader_timstof_recovered_isos <- bader_timstof[
        match(recovered_prots_isoforms, rownames(bader_timstof)),
]

bader_timstof_recovered_isos_dups <- gsub(
        "\\-.*", "", rownames(bader_timstof_recovered_isos)
)[
        duplicated(gsub("\\-.*", "", rownames(bader_timstof_recovered_isos)))
]

recovered_prots_isoforms_dups <- recovered_prots_isoforms[
        gsub("\\-.*", "", recovered_prots_isoforms) %in% bader_timstof_recovered_isos_dups
]

bader_timstof_recovered_isos_dups_df <- bader_timstof[
        match(recovered_prots_isoforms_dups, rownames(bader_timstof)),
]

bader_timstof_recovered_isos_long <- bader_timstof_recovered_isos_dups_df %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()


ggplot(bader_timstof_recovered_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()


bader_timstof_no_dups <- as.matrix(bader_timstof)

rownames(bader_timstof_no_dups)[
        rownames(bader_timstof_no_dups) %in% recovered_prots_isoforms
] <- gsub("\\-.*", "", rownames(bader_timstof_no_dups)[
        rownames(bader_timstof_no_dups) %in% recovered_prots_isoforms
])

rownames(bader_timstof_no_dups) <- gsub("\\-.*", "", rownames(bader_timstof_no_dups))

bader_timstof_no_dups <- t(
        sapply(
                unique(rownames(bader_timstof_no_dups)),
                function(x) colMeans(
                        bader_timstof_no_dups[x, , drop = FALSE], na.rm = TRUE
                )
        )
)

bader_timstof_no_dups[is.nan(bader_timstof_no_dups)] <- NA

recovered_prots_isoforms_dups


yet_missing <- mapping$Accession[
        !mapping$Accession %in% rownames(bader_timstof_no_dups)
]


bader_mapping_mart_f <- file.path(out_dir, "bader_mapping.csv")
if (!file.exists(bader_mapping_mart_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        bader_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = rownames(bader_timstof_no_dups),
                mart = mart
        )
        write.csv(bader_mapping, file = bader_mapping_mart_f)
}else{
        bader_mapping <- read.csv(bader_mapping_mart_f)
}

map[map$uniprotswissprot %in% yet_missing, ]$hgnc_symbol %in% bader_mapping$hgnc_symbol


################################################################################
# EMORY 3
################################################################################

emory_3_raw_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_3/Master.3platform.2fluid/1b.PD.TMT_rawAbundances.NoCorrection.csv"


emory_3_raw <- read.csv(emory_3_raw_f, row.names = 1)

plot(density(as.matrix(emory_3_raw), na.rm = TRUE))

dim(emory_3_raw)
emory_3_raw[1:10, 1:10]


aligned_out <- load("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/Emory_3/Master.3platform.2fluid/2.Aligned.outputs.RData")
dim(emory_3_raw)

emory_3_csf_prot <- cleanDat.MS.CSFneat

plot(density(as.matrix(emory_3_csf_prot), na.rm = TRUE))

dim(as.matrix(cleanDat.MS.CSFdepl))

dim(rawRFUdat.Soma.CSFneat.36)
plot(density(rawRFUdat.Soma.CSFneat.36, na.rm = TRUE))

emory_3_csf_mdat <- numericMeta.MS

emory_3_uniprot_IDs <- sapply(rownames(emory_3_csf_prot), function(x) strsplit(x, split = "|", fixed = TRUE)[[1]][2])

length(intersect(emory_3_uniprot_IDs, mapping$Accession))/length(emory_3_uniprot_IDs)

length(emory_3_uniprot_IDs[!emory_3_uniprot_IDs %in% mapping$Accession])

emory_3_uniprot_IDs_missing <- emory_3_uniprot_IDs[!emory_3_uniprot_IDs %in% mapping$Accession]

length(emory_3_uniprot_IDs %>% intersect(mapping$Accession) %>% intersect(rownames(bader_timstof_no_dups)))


emory_3_mapping_f <- file.path(out_dir, "emory_3_mapping.csv")
if (!file.exists(emory_3_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        emory_3_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = unique(gsub("\\-.*", "", emory_3_uniprot_IDs)),
                mart = mart
        )
        write.csv(emory_3_mapping, file = emory_3_mapping_f)
}else{
        emory_3_mapping <- read.csv(emory_3_mapping_f)
}


emory_3_mapping[emory_3_mapping$uniprotswissprot %in% emory_3_uniprot_IDs_missing, ]$hgnc_symbol %in% bader_mapping$hgnc_symbol
sum(emory_3_mapping[emory_3_mapping$uniprotswissprot %in% emory_3_uniprot_IDs_missing, ]$hgnc_symbol %in% map$hgnc_symbol)

emory_3_uniprotIDs <- sapply(
        rownames(emory_3_csf_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
emory_3_uniprotIDs_deisoformed <- gsub("\\-.*", "", emory_3_uniprotIDs)
emory_3_uniprotIDs_deisoformed_dups <- emory_3_uniprotIDs_deisoformed[duplicated(emory_3_uniprotIDs_deisoformed)]
emory_3_uniprotIDs_deisoformed <- unique(emory_3_uniprotIDs_deisoformed)


emory_3_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(emory_3_csf_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% emory_3_uniprotIDs_deisoformed_dups


emory_3_csf_prot_dup_isos_long <- emory_3_csf_prot[emory_3_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(emory_3_csf_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform, if it's present. If it's not, we compute the mean.
emory_3_csf_prot_no_isos <- lapply(emory_3_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(emory_3_csf_prot))) == 1) {
                emory_3_csf_prot[grepl(x, rownames(emory_3_csf_prot)), ]
        } else if (any(sapply(rownames(emory_3_csf_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x)){
                emory_3_csf_prot[
                        sapply(rownames(emory_3_csf_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x, 
                ]
        } else {
                colMeans(emory_3_csf_prot[grepl(x, rownames(emory_3_csf_prot)), ], na.rm = TRUE)
        }
})

emory_3_csf_prot_no_isos <- do.call(rbind, emory_3_csf_prot_no_isos)
rownames(emory_3_csf_prot_no_isos) <- emory_3_uniprotIDs_deisoformed

################################################################################
# syn20821165
################################################################################

syn20821165_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/syn20821165"

list.files(file.path(syn20821165_dir, "CSF_discovery_cohort_1_asymAD"))

syn20821165_coh_1_prot <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_1_asymAD/CSF Replication Cohort 1-cleanDat.xlsx"
                )
        )
)
dim(syn20821165_coh_1_prot)

rownames(syn20821165_coh_1_prot) <- syn20821165_coh_1_prot$Symbol.UniprotID

syn20821165_coh_1_prot_names <- rownames(syn20821165_coh_1_prot)

syn20821165_coh_1_prot <- syn20821165_coh_1_prot[, 2:ncol(syn20821165_coh_1_prot)]

syn20821165_coh_1_prot <- apply(
        syn20821165_coh_1_prot, 2, function(x) ifelse(x == "NA", NA, x)
)

syn20821165_coh_1_prot <- apply(syn20821165_coh_1_prot, 2, as.numeric)
rownames(syn20821165_coh_1_prot) <- syn20821165_coh_1_prot_names


plot(density(as.matrix(syn20821165_coh_1_prot), na.rm = TRUE))

syn20821165_coh_1_prot_info <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_1_asymAD/CSF Replication Cohort 1-search result.xlsx"
                )
        )
)

syn20821165_coh_1_mdat <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_1_asymAD/CSF Replication Cohort 1-sample traits.xlsx"
                )
        )
)


syn20821165_coh_1_uniprotIDs <- sapply(
        rownames(syn20821165_coh_1_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20821165_coh_1_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20821165_coh_1_uniprotIDs)
syn20821165_coh_1_uniprotIDs_deisoformed_dups <- syn20821165_coh_1_uniprotIDs_deisoformed[duplicated(syn20821165_coh_1_uniprotIDs_deisoformed)]
syn20821165_coh_1_uniprotIDs_deisoformed <- unique(syn20821165_coh_1_uniprotIDs_deisoformed)


syn20821165_coh_1_mapping_f <- file.path(out_dir, "syn20821165_coh_1_mapping.csv")
if (!file.exists(syn20821165_coh_1_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20821165_coh_1_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20821165_coh_1_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20821165_coh_1_mapping, file = syn20821165_coh_1_mapping_f)
}else{
        syn20821165_coh_1_mapping <- read.csv(syn20821165_coh_1_mapping_f)
}

syn20821165_coh_1_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20821165_coh_1_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20821165_coh_1_uniprotIDs_deisoformed_dups
        

syn20821165_coh_1_prot_dup_isos_long <- syn20821165_coh_1_prot[syn20821165_coh_1_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20821165_coh_1_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

list.files(file.path(syn20821165_dir, "CSF_discovery_cohort_2_FNIH"))


# Here it looks like the isoforms are different. Just keep the one with no
# suffix

sapply(rownames(syn20821165_coh_1_prot), function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2])

sum(grepl('-', syn20821165_coh_1_uniprotIDs))

# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform.
syn20821165_coh_1_prot_no_isos <- lapply(syn20821165_coh_1_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20821165_coh_1_prot))) == 1) {
                syn20821165_coh_1_prot[grepl(x, rownames(syn20821165_coh_1_prot)), ]
        } else {
                syn20821165_coh_1_prot[
                        sapply(rownames(syn20821165_coh_1_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x, 
                ]
        }
})

syn20821165_coh_1_prot_no_isos <- do.call(rbind, syn20821165_coh_1_prot_no_isos)
rownames(syn20821165_coh_1_prot_no_isos) <- syn20821165_coh_1_uniprotIDs_deisoformed









#syn20821165_coh_2_prot <- data.frame(
#        read.csv(
#                file.path(
#                        syn20821165_dir,
#                        "CSF_discovery_cohort_2_FNIH/2b.Unregressed_Batch-corrected_cleanDat.csv"
#                ),
#                row.names = 1
#        )
#)

syn20821165_coh_2_prot <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        #"CSF_discovery_cohort_2_FNIH/2.Unregressed_Batch-corrected_PD-normalized_TMT_reporter_abundance.csv" # This is un-ransformed
                        "CSF_discovery_cohort_2_FNIH/2b.Unregressed_Batch-corrected_cleanDat.csv" # This one looks as age info has not been regressed out (SD of cor with age is ~0.1)
                        #"CSF_discovery_cohort_2_FNIH/3.cleanDat.csv" # This one looks as age info has been regressed out (SD of cor with age is ~0.02)
                ),
                row.names = 1
        )
)
dim(syn20821165_coh_2_prot)

plot(density(as.matrix(syn20821165_coh_2_prot), na.rm = TRUE))

syn20821165_coh_2_prot_info <- data.frame(
        read.table(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_2_FNIH/1.PD-proteins.txt"
                ),
                sep = '\t',
                header = TRUE
        )
)

syn20821165_coh_2_mdat <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_2_FNIH/0.Traits.csv"
                ),
                row.names = 1
        )
)

syn20821165_coh_2_uniprotIDs <- sapply(
        rownames(syn20821165_coh_2_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20821165_coh_2_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20821165_coh_2_uniprotIDs)
syn20821165_coh_2_uniprotIDs_deisoformed_dups <- syn20821165_coh_2_uniprotIDs_deisoformed[duplicated(syn20821165_coh_2_uniprotIDs_deisoformed)]
syn20821165_coh_2_uniprotIDs_deisoformed <- unique(syn20821165_coh_2_uniprotIDs_deisoformed)


syn20821165_coh_2_mapping_f <- file.path(out_dir, "syn20821165_coh_2_mapping.csv")
if (!file.exists(syn20821165_coh_2_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20821165_coh_2_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20821165_coh_2_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20821165_coh_2_mapping, file = syn20821165_coh_2_mapping_f)
}else{
        syn20821165_coh_2_mapping <- read.csv(syn20821165_coh_2_mapping_f)
}

syn20821165_coh_2_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20821165_coh_2_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20821165_coh_2_uniprotIDs_deisoformed_dups


syn20821165_coh_2_prot_dup_isos_long <- syn20821165_coh_2_prot[syn20821165_coh_2_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20821165_coh_2_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

list.files(file.path(syn20821165_dir, "CSF_discovery_cohort_2_FNIH"))


# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform. For P10636 (tau protein), keep -MTBR version
syn20821165_coh_2_prot_no_isos <- lapply(syn20821165_coh_2_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20821165_coh_2_prot))) == 1) {
                syn20821165_coh_2_prot[grepl(x, rownames(syn20821165_coh_2_prot)), ]
        } else {
                syn20821165_coh_2_prot[
                        sapply(rownames(syn20821165_coh_2_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x | sapply(rownames(syn20821165_coh_2_prot), function(r) strsplit(r, "\\|")[[1]][2]) == paste0(x, "-MTBR"), 
                ]
        }
})

syn20821165_coh_2_prot_no_isos <- do.call(rbind, syn20821165_coh_2_prot_no_isos)
rownames(syn20821165_coh_2_prot_no_isos) <- syn20821165_coh_2_uniprotIDs_deisoformed











syn20821165_coh_3_prot <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_3_Nestle_Dayon/cleanDat-batchCorrected_Averaged_Technical_Replicates.csv"
                ),
                row.names = 1
        )
)
dim(syn20821165_coh_3_prot)

plot(density(as.matrix(syn20821165_coh_3_prot), na.rm = TRUE))

syn20821165_coh_3_prot_info <- data.frame(
        read.table(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_3_Nestle_Dayon/PD2_DayonCSF_Proteins.txt"
                ),
                sep = '\t',
                header = TRUE
        )
)

syn20821165_coh_3_mdat <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_3_Nestle_Dayon/DayonCSFreplicates_singleShotTMT_traits.csv"
                ),
                row.names = 1
        )
)


syn20821165_coh_3_uniprotIDs <- sapply(
        rownames(syn20821165_coh_3_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20821165_coh_3_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20821165_coh_3_uniprotIDs)
syn20821165_coh_3_uniprotIDs_deisoformed_dups <- syn20821165_coh_3_uniprotIDs_deisoformed[duplicated(syn20821165_coh_3_uniprotIDs_deisoformed)]
syn20821165_coh_3_uniprotIDs_deisoformed <- unique(syn20821165_coh_3_uniprotIDs_deisoformed)

syn20821165_coh_3_mapping_f <- file.path(out_dir, "syn20821165_coh_3_mapping.csv")
if (!file.exists(syn20821165_coh_3_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20821165_coh_3_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20821165_coh_3_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20821165_coh_3_mapping, file = syn20821165_coh_3_mapping_f)
}else{
        syn20821165_coh_3_mapping <- read.csv(syn20821165_coh_3_mapping_f)
}

syn20821165_coh_3_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20821165_coh_3_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20821165_coh_3_uniprotIDs_deisoformed_dups


syn20821165_coh_3_prot_dup_isos_long <- syn20821165_coh_3_prot[syn20821165_coh_3_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20821165_coh_3_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform. For P10636 (tau protein), keep -MTBR version
syn20821165_coh_3_prot_no_isos <- lapply(syn20821165_coh_3_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20821165_coh_3_prot))) == 1) {
                syn20821165_coh_3_prot[grepl(x, rownames(syn20821165_coh_3_prot)), ]
        } else {
                syn20821165_coh_3_prot[
                        sapply(rownames(syn20821165_coh_3_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x, 
                ]
        }
})

syn20821165_coh_3_prot_no_isos <- do.call(rbind, syn20821165_coh_3_prot_no_isos)
rownames(syn20821165_coh_3_prot_no_isos) <- syn20821165_coh_3_uniprotIDs_deisoformed




syn20821165_coh_4_prot <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_4_Multi_disease/CSF Replication Cohort 4-cleanDat.csv"
                ),
                row.names = 1
        )
)
dim(syn20821165_coh_4_prot)

plot(density(as.matrix(syn20821165_coh_4_prot), na.rm = TRUE))

syn20821165_coh_4_prot_info <- data.frame(
        read.table(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_4_Multi_disease/CSF Replication Cohort 4-PD_round1multiconcensus_nosnfilter_normabundance_Proteins.txt"
                ),
                sep = '\t',
                header = TRUE
        )
)

syn20821165_coh_4_mdat <- data.frame(
        read.csv(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_4_Multi_disease/CSF Replication Cohort 4-traits.csv"
                ),
                row.names = 1
        )
)


syn20821165_coh_4_uniprotIDs <- sapply(
        rownames(syn20821165_coh_4_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20821165_coh_4_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20821165_coh_4_uniprotIDs)
syn20821165_coh_4_uniprotIDs_deisoformed_dups <- syn20821165_coh_4_uniprotIDs_deisoformed[duplicated(syn20821165_coh_4_uniprotIDs_deisoformed)]
syn20821165_coh_4_uniprotIDs_deisoformed <- unique(syn20821165_coh_4_uniprotIDs_deisoformed)


syn20821165_coh_4_mapping_f <- file.path(out_dir, "syn20821165_coh_4_mapping.csv")
if (!file.exists(syn20821165_coh_4_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20821165_coh_4_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20821165_coh_4_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20821165_coh_4_mapping, file = syn20821165_coh_4_mapping_f)
}else{
        syn20821165_coh_4_mapping <- read.csv(syn20821165_coh_4_mapping_f)
}


syn20821165_coh_4_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20821165_coh_4_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20821165_coh_4_uniprotIDs_deisoformed_dups


syn20821165_coh_4_prot_dup_isos_long <- syn20821165_coh_4_prot[syn20821165_coh_4_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20821165_coh_4_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform.
syn20821165_coh_4_prot_no_isos <- lapply(syn20821165_coh_4_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20821165_coh_4_prot))) == 1) {
                syn20821165_coh_4_prot[grepl(x, rownames(syn20821165_coh_4_prot)), ]
        } else {
                syn20821165_coh_4_prot[
                        sapply(rownames(syn20821165_coh_4_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x, 
                ]
        }
})

syn20821165_coh_4_prot_no_isos <- do.call(rbind, syn20821165_coh_4_prot_no_isos)
rownames(syn20821165_coh_4_prot_no_isos) <- syn20821165_coh_4_uniprotIDs_deisoformed






syn20821165_csf_disc_coh_fracTMT_prot <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_fractionated_TMT/clean data_discovery AD CSF 40_fractionated TMT.xlsx"
                )
        )
)
dim(syn20821165_csf_disc_coh_fracTMT_prot)


rownames(syn20821165_csf_disc_coh_fracTMT_prot) <- syn20821165_csf_disc_coh_fracTMT_prot$Symbol.UniprotID

syn20821165_csf_disc_coh_fracTMT_prot_names <- rownames(syn20821165_csf_disc_coh_fracTMT_prot)

syn20821165_csf_disc_coh_fracTMT_prot <- syn20821165_csf_disc_coh_fracTMT_prot[, 2:ncol(syn20821165_csf_disc_coh_fracTMT_prot)]

syn20821165_csf_disc_coh_fracTMT_prot <- apply(
        syn20821165_csf_disc_coh_fracTMT_prot, 2, function(x) ifelse(x == "NA", NA, x)
)

syn20821165_csf_disc_coh_fracTMT_prot <- apply(syn20821165_csf_disc_coh_fracTMT_prot, 2, as.numeric)
rownames(syn20821165_csf_disc_coh_fracTMT_prot) <- syn20821165_csf_disc_coh_fracTMT_prot_names

plot(density(syn20821165_csf_disc_coh_fracTMT_prot, na.rm = TRUE))

syn20821165_csf_disc_coh_fracTMT_prot_info <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_fractionated_TMT/searching result_discovery AD CSF 40_fractionated TMT.xlsx"
                )
        )
)

syn20821165_csf_disc_coh_fracTMT_mdat <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20821165_dir,
                        "CSF_discovery_cohort_fractionated_TMT/sample trait_discovery AD CSF 40_fractionated TMT.xlsx"
                )
        )
)

syn20821165_csf_disc_coh_fracTMT_uniprotIDs <- sapply(
        rownames(syn20821165_csf_disc_coh_fracTMT_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20821165_csf_disc_coh_fracTMT_uniprotIDs)
syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed_dups <- syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed[duplicated(syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed)]
syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed <- unique(syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed)


syn20821165_csf_disc_coh_fracTMT_mapping_f <- file.path(out_dir, "syn20821165_csf_disc_coh_fracTMT_mapping.csv")
if (!file.exists(syn20821165_csf_disc_coh_fracTMT_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20821165_csf_disc_coh_fracTMT_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20821165_csf_disc_coh_fracTMT_mapping, file = syn20821165_csf_disc_coh_fracTMT_mapping_f)
}else{
        syn20821165_csf_disc_coh_fracTMT_mapping <- read.csv(syn20821165_csf_disc_coh_fracTMT_mapping_f)
}

syn20821165_csf_disc_coh_fracTMT_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20821165_csf_disc_coh_fracTMT_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed_dups


syn20821165_csf_disc_coh_fracTMT_prot_dup_isos_long <- syn20821165_csf_disc_coh_fracTMT_prot[syn20821165_csf_disc_coh_fracTMT_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20821165_csf_disc_coh_fracTMT_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. When there are multiple, we keep the canonical
# isoform, if it's present. If it's not, we compute the mean.
syn20821165_csf_disc_coh_fracTMT_prot_no_isos <- lapply(syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20821165_csf_disc_coh_fracTMT_prot))) == 1) {
                syn20821165_csf_disc_coh_fracTMT_prot[grepl(x, rownames(syn20821165_csf_disc_coh_fracTMT_prot)), ]
        } else if (any(sapply(rownames(syn20821165_csf_disc_coh_fracTMT_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x)){
                syn20821165_csf_disc_coh_fracTMT_prot[
                        sapply(rownames(syn20821165_csf_disc_coh_fracTMT_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x, 
                ]
        } else {
                colMeans(syn20821165_csf_disc_coh_fracTMT_prot[grepl(x, rownames(syn20821165_csf_disc_coh_fracTMT_prot)), ], na.rm = TRUE)
        }
})

syn20821165_csf_disc_coh_fracTMT_prot_no_isos <- do.call(rbind, syn20821165_csf_disc_coh_fracTMT_prot_no_isos)
rownames(syn20821165_csf_disc_coh_fracTMT_prot_no_isos) <- syn20821165_csf_disc_coh_fracTMT_uniprotIDs_deisoformed


################################################################################
# syn20933797
################################################################################

syn20933797_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/syn20933797"

syn20933797_coh_1_prot <- data.frame(
        read.csv(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_1_sample_FNIH_sympt_AD_and_ctrls/2.Unregressed_Batch-corrected_PD-normalized_TMT_reporter_abundance.csv"
                        #"CSF_cohort_1_sample_FNIH_sympt_AD_and_ctrls/3.cleanDat.csv"
                        #"CSF_cohort_1_sample_FNIH_sympt_AD_and_ctrls/2b.Unregressed_Batch-corrected_cleanDat.csv"
                ),
                row.names = 1
        )
)
dim(syn20933797_coh_1_prot)

plot(density(as.matrix(syn20933797_coh_1_prot), na.rm = TRUE))

syn20933797_coh_1_prot_info <- data.frame(
        read.table(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_1_sample_FNIH_sympt_AD_and_ctrls/1.PD-proteins.txt"
                ),
                sep = '\t',
                header = TRUE
        )
)

syn20933797_coh_1_mdat <- data.frame(
        read.csv(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_1_sample_FNIH_sympt_AD_and_ctrls/0.Traits.csv"
                ),
                row.names = 1
        )
)


syn20933797_coh_1_uniprotIDs <- sapply(
        rownames(syn20933797_coh_1_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20933797_coh_1_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20933797_coh_1_uniprotIDs)
syn20933797_coh_1_uniprotIDs_deisoformed_dups <- syn20933797_coh_1_uniprotIDs_deisoformed[duplicated(syn20933797_coh_1_uniprotIDs_deisoformed)]
syn20933797_coh_1_uniprotIDs_deisoformed <- unique(syn20933797_coh_1_uniprotIDs_deisoformed)


syn20933797_coh_1_mapping_f <- file.path(out_dir, "syn20933797_coh_1_mapping.csv")
if (!file.exists(syn20933797_coh_1_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20933797_coh_1_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20933797_coh_1_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20933797_coh_1_mapping, file = syn20933797_coh_1_mapping_f)
}else{
        syn20933797_coh_1_mapping <- read.csv(syn20933797_coh_1_mapping_f)
}


syn20933797_coh_1_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20933797_coh_1_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20933797_coh_1_uniprotIDs_deisoformed_dups


syn20933797_coh_1_prot_dup_isos_long <- syn20933797_coh_1_prot[syn20933797_coh_1_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20933797_coh_1_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. For P10636 (tau), keep the MTBR one
syn20933797_coh_1_prot_no_isos <- lapply(syn20933797_coh_1_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20933797_coh_1_prot))) == 1) {
                syn20933797_coh_1_prot[grepl(x, rownames(syn20933797_coh_1_prot)), ]
        } else {
                syn20933797_coh_1_prot[
                        sapply(rownames(syn20933797_coh_1_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x | sapply(rownames(syn20933797_coh_1_prot), function(r) strsplit(r, "\\|")[[1]][2]) == paste0(x, "-MTBR"), 
                ]
        }
})

syn20933797_coh_1_prot_no_isos <- do.call(rbind, syn20933797_coh_1_prot_no_isos)
rownames(syn20933797_coh_1_prot_no_isos) <- syn20933797_coh_1_uniprotIDs_deisoformed











syn20933797_coh_2_prot <- data.frame(
        read.csv(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_2_sample_biomarker_defined_asympt_AD_sympt_AD_and_ctrls/2.Unregressed_Batch-corrected_PD-TMT_reporter_abundance.csv"
                        #"CSF_cohort_2_sample_biomarker_defined_asympt_AD_sympt_AD_and_ctrls/2b.unregressed_Batch-corrected_cleanDat.csv"
                        #"CSF_cohort_2_sample_biomarker_defined_asympt_AD_sympt_AD_and_ctrls/3.cleanDat.csv"
                ),
                row.names = 1
        )
)
dim(syn20933797_coh_2_prot)

plot(density(as.matrix(syn20933797_coh_2_prot), na.rm = TRUE))

syn20933797_coh_2_prot_info <- data.frame(
        readxl::read_xlsx(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_2_sample_biomarker_defined_asympt_AD_sympt_AD_and_ctrls/1.PD-proteins.txt.xlsx"
                )
        )
)

syn20933797_coh_2_mdat <- data.frame(
        read.csv(
                file.path(
                        syn20933797_dir,
                        "CSF_cohort_2_sample_biomarker_defined_asympt_AD_sympt_AD_and_ctrls/0.Traits.csv"
                ),
                row.names = 1
        )
)


syn20933797_coh_2_uniprotIDs <- sapply(
        rownames(syn20933797_coh_2_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)
syn20933797_coh_2_uniprotIDs_deisoformed <- gsub("\\-.*", "", syn20933797_coh_2_uniprotIDs)
syn20933797_coh_2_uniprotIDs_deisoformed_dups <- syn20933797_coh_2_uniprotIDs_deisoformed[duplicated(syn20933797_coh_2_uniprotIDs_deisoformed)]
syn20933797_coh_2_uniprotIDs_deisoformed <- unique(syn20933797_coh_2_uniprotIDs_deisoformed)


syn20933797_coh_2_mapping_f <- file.path(out_dir, "syn20933797_coh_2_mapping.csv")
if (!file.exists(syn20933797_coh_2_mapping_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        syn20933797_coh_2_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = syn20933797_coh_2_uniprotIDs_deisoformed,
                mart = mart
        )
        write.csv(syn20933797_coh_2_mapping, file = syn20933797_coh_2_mapping_f)
}else{
        syn20933797_coh_2_mapping <- read.csv(syn20933797_coh_2_mapping_f)
}


syn20933797_coh_2_is_iso_vec <- gsub("\\-.*", "", sapply(
        rownames(syn20933797_coh_2_prot),
        function(x) strsplit(x, split = '|', fixed = TRUE)[[1]][2]
)) %in% syn20933797_coh_2_uniprotIDs_deisoformed_dups


syn20933797_coh_2_prot_dup_isos_long <- syn20933797_coh_2_prot[syn20933797_coh_2_is_iso_vec, ] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()

ggplot(syn20933797_coh_2_prot_dup_isos_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

# Here it looks like the isoforms are different. Just keep the one with no
# suffix


# Here, we iterate through each unique protein ID and, where there's only
# one match, we keep that vector. For P10636 (tau), keep the MTBR one
syn20933797_coh_2_prot_no_isos <- lapply(syn20933797_coh_2_uniprotIDs_deisoformed, function(x) {
        if (sum(grepl(x, rownames(syn20933797_coh_2_prot))) == 1) {
                syn20933797_coh_2_prot[grepl(x, rownames(syn20933797_coh_2_prot)), ]
        } else {
                syn20933797_coh_2_prot[
                        sapply(rownames(syn20933797_coh_2_prot), function(r) strsplit(r, "\\|")[[1]][2]) == x | sapply(rownames(syn20933797_coh_2_prot), function(r) strsplit(r, "\\|")[[1]][2]) == paste0(x, "-MTBR"), 
                ]
        }
})

syn20933797_coh_2_prot_no_isos <- do.call(rbind, syn20933797_coh_2_prot_no_isos)
rownames(syn20933797_coh_2_prot_no_isos) <- syn20933797_coh_2_uniprotIDs_deisoformed

################################################################################
# PRIDE-AD
################################################################################

pride_ad_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/PRIDE_AD/Discovery_cohort/DataSet_PRIDE_ADstudy.csv"
pride_ad <- read.csv(pride_ad_f, row.names = 1)

pride_ad_mdat <- pride_ad[, 1:3]

pride_ad_mdat$individualID <- rownames(pride_ad_mdat)

pride_ad_prot <- pride_ad[, 4:ncol(pride_ad)]

colnames(pride_ad_prot)

olink_env <- as.list(getNamespace("OlinkAnalyze"))

panel_maps <- olink_env$olink_wide_spec

library(rvest)

# Get prtoein-Uniprot equivalence matrices from olinks websites
card_II_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-cardiovascular-ii-panel.html"
card_III_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-cardiovascular-iii-panel.html"
cardimet_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-cardiometabolic-panel.html"
cellregu_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-cell-regulation-panel.html"
developm_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-development-panel.html"
immunres_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-immune-response-panel.html"
inflammt_url <- "https://olinkpanel.creative-proteomics.com/table-list-of-olink-explore-384-inflammation-panel.html"
metabolm_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-metabolism-panel.html"
neurolog_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-neurology-panel.html"
oncoloII_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-oncology-ii-panel.html"
organdam_url <- "https://olinkpanel.creative-proteomics.com/panel/olink-target-96-organ-damage-panel.html"


card_II_tabs <- read_html(card_II_url) |> html_table(fill = TRUE)
card_II_tabs <- as.data.frame(card_II_tabs[[1]])
colnames(card_II_tabs) <- unlist(as.vector(card_II_tabs[1, ]))
card_II_tabs <- card_II_tabs[2:nrow(card_II_tabs), ]
card_II_tabs$panel <- "Cardiovascular II"

card_III_tabs <- read_html(card_III_url) |> html_table(fill = TRUE)
card_III_tabs <- as.data.frame(card_III_tabs[[1]])
colnames(card_III_tabs) <- unlist(as.vector(card_III_tabs[1, ]))
card_III_tabs <- card_III_tabs[2:nrow(card_III_tabs), ]
card_III_tabs$panel <- "Cardiovascular III"

cardimet_tabs <- read_html(cardimet_url) |> html_table(fill = TRUE)
cardimet_tabs <- as.data.frame(cardimet_tabs[[1]])
colnames(cardimet_tabs) <- unlist(as.vector(cardimet_tabs[1, ]))
cardimet_tabs <- cardimet_tabs[2:nrow(cardimet_tabs), ]
cardimet_tabs$panel <- "Cardiometabolic"

cellregu_tabs <- read_html(cellregu_url) |> html_table(fill = TRUE)
cellregu_tabs <- as.data.frame(cellregu_tabs[[1]])
colnames(cellregu_tabs) <- unlist(as.vector(cellregu_tabs[1, ]))
cellregu_tabs <- cellregu_tabs[2:nrow(cellregu_tabs), ]
cellregu_tabs$panel <- "Cell regulation"

developm_tabs <- read_html(developm_url) |> html_table(fill = TRUE)
developm_tabs <- as.data.frame(developm_tabs[[1]])
colnames(developm_tabs) <- unlist(as.vector(developm_tabs[1, ]))
developm_tabs <- developm_tabs[2:nrow(developm_tabs), ]
developm_tabs$panel <- "Development"
colnames(developm_tabs) <- c(
        "Protein Category", "UniProt ID", "Gene", "Protein Name", "panel"
)

immunres_tabs <- read_html(immunres_url) |> html_table(fill = TRUE)
immunres_tabs <- as.data.frame(immunres_tabs[[1]])
colnames(immunres_tabs) <- unlist(as.vector(immunres_tabs[1, ]))
immunres_tabs <- immunres_tabs[2:nrow(immunres_tabs), ]
immunres_tabs$panel <- "Immune response"

inflammt_tabs <- read_html(inflammt_url) |> html_table(fill = TRUE)
inflammt_tabs <- as.data.frame(inflammt_tabs[[1]])
colnames(inflammt_tabs) <- unlist(as.vector(inflammt_tabs[1, ]))
inflammt_tabs <- inflammt_tabs[2:nrow(inflammt_tabs), ]
inflammt_tabs$panel <- "Inflammation"
colnames(inflammt_tabs) <- gsub("Protein name", "Protein Name", colnames(inflammt_tabs))
inflammt_tabs$`Protein Category` <- NA
inflammt_tabs <- inflammt_tabs[, colnames(immunres_tabs)]

metabolm_tabs <- read_html(metabolm_url) |> html_table(fill = TRUE)
metabolm_tabs <- as.data.frame(metabolm_tabs[[1]])
colnames(metabolm_tabs) <- unlist(as.vector(metabolm_tabs[1, ]))
metabolm_tabs <- metabolm_tabs[2:nrow(metabolm_tabs), ]
metabolm_tabs$panel <- "Metabolism"

neurolog_tabs <- read_html(neurolog_url) |> html_table(fill = TRUE)
neurolog_tabs <- as.data.frame(neurolog_tabs[[1]])
colnames(neurolog_tabs) <- unlist(as.vector(neurolog_tabs[1, ]))
neurolog_tabs <- neurolog_tabs[2:nrow(neurolog_tabs), ]
neurolog_tabs$panel <- "Neurology"

oncoloII_tabs <- read_html(oncoloII_url) |> html_table(fill = TRUE)
oncoloII_tabs <- as.data.frame(oncoloII_tabs[[1]])
colnames(oncoloII_tabs) <- unlist(as.vector(oncoloII_tabs[1, ]))
oncoloII_tabs <- oncoloII_tabs[2:nrow(oncoloII_tabs), ]
oncoloII_tabs$panel <- "Oncology II"

organdam_tabs <- read_html(organdam_url) |> html_table(fill = TRUE)
organdam_tabs <- as.data.frame(organdam_tabs[[1]])
colnames(organdam_tabs) <- unlist(as.vector(organdam_tabs[1, ]))
organdam_tabs <- organdam_tabs[2:nrow(organdam_tabs), ]
organdam_tabs$panel <- "Organ damage"

olink_panel_mapping <- rbind.data.frame(
        card_II_tabs,
        card_III_tabs,
        cardimet_tabs,
        cellregu_tabs,
        developm_tabs,
        immunres_tabs,
        inflammt_tabs,
        metabolm_tabs,
        neurolog_tabs,
        oncoloII_tabs,
        organdam_tabs
)
write.csv(olink_panel_mapping, file = file.path(out_dir, "pride_ad_olink_mapping.csv"))

colnames(pride_ad_prot)[!colnames(pride_ad_prot) %in% make.names(olink_panel_mapping$Gene)]

if (!require(pdftools)) install.packages("pdftools")
library(pdftools)
library(stringr)

txt <- pdf_text("/Users/guillem.santamaria/Downloads/NIHMS1902789-supplement-SupplementaryTables1-4.pdf")
lines <- unlist(strsplit(txt, "\n"))[3:654]
lines_splits <- sapply(lines, function(x) strsplit(x, split = "\\s{2,}"))

lines_splits[which(as.vector(lapply(lines_splits, length) == 5))]

# Do dataframe of the ones that have the number of items they should have,
# and then explore the ones that were not parsed correctly.
mapping_olink <- do.call(rbind.data.frame, lines_splits[which(as.vector(lapply(lines_splits, length) == 5))])
colnames(mapping_olink) <- c("Name", "Abbreviation", "Uniprot", "Panel", "CSF_detectability")
mapping_olink <- mapping_olink[2:nrow(mapping_olink), ]

lines[which(as.vector(lapply(lines_splits, length) != 5))][3]

unmapped_manual <- data.frame(
        Name = c(
                "Tumor necrosis factor receptor superfamily member 10A",
                "Vascular endothelial growth factor D",
                "Interleukin-6",
                "Granulocyte-macrophage colony-stimulating factor receptor subunit alpha",
                "Stem cell factor",
                "Interleukin-12",
                "CMP-N-acetylneuraminate-beta-galactosamide-alpha-2,3-sialyltransferase 1",
                "MHC class I polypeptide-related sequence A/B",
                "Interleukin-27",
                "MAM domain-containing glycosylphosphatidylinositol anchor protein 1",
                "WAP, Kazal, immunoglobulin, Kunitz and NTR domain-containing protein 2",
                "WAP, Kazal, immunoglobulin, Kunitz and NTR domain-containing protein 1",
                "Phosphoprotein associated with glycosphingolipid-enriched microdomains 1",
                "Galactoside 3(4)-L-fucosyltransferase,Alpha-(1,3)-fucosyltransferase 3/5"
        ),
        Abbreviation = c(
                "TNFRSF10A",
                "VEGFD",
                "IL6",
                "GM-CSF-R-alpha",
                "SCF",
                "IL-12",
                "ST3GAL1",
                "MIC-A/B",
                "IL-27",
                "MDGA1",
                "WFIKKN2",
                "WFIKKN1",
                "PAG1",
                "FUT3/FUT5"
        ),
        Uniprot = c(
                "O00220",
                "O43915",
                "P05231",
                "P15509",
                "P21583",
                "P29460,P29459",
                "Q11201",
                "Q29983,Q29980",
                "Q8NEV9,Q14213",
                "Q8NFP4",
                "Q8TEU8",
                "Q96NZ8",
                "Q9NWQ8",
                "Q11128,P21217"
        ),
        Panel = c(
                "CVD II, Cell-Reg",
                "CVD II, Cell-Reg",
                "Inf, Onc II, CVD II, IR",
                "Neuro",
                "Inf, CVD II, Onc II",
                "Neuro",
                "OD",
                "Onc II",
                "CVD II",
                "Neuro",
                "Dev",
                "Neuro",
                "Met",
                "Dev"
        ),
        CSF_detectability = c(
                "100%; 100%",
                "100%; 99.79%",
                "100%; 100%; 100%; 100%",
                "100%",
                "100%; 100%",
                "99.90%",
                "99.69%",
                "99.59%",
                "100%",
                "100%",
                "100%",
                "92.25%",
                "100%",
                "95.86%"
        )
)


mapping_olink <- rbind.data.frame(mapping_olink, unmapped_manual)

#mapping_olink_deepseek <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/PRIDE_AD/olink_mapping_deepseek.csv")

#mapping_olink_deepseek[!mapping_olink_deepseek$Matrix_Column_Name %in% make.names(mapping_olink$Abbreviation), ]

mapping_olink[grep("IL", mapping_olink$Abbreviation), ]
mapping_olink[grep("TGF", mapping_olink$Abbreviation), ]
colnames(pride_ad_prot)[grep("IL", colnames(pride_ad_prot))]


colnames(pride_ad_prot)[!colnames(pride_ad_prot) %in% make.names(mapping_olink$Abbreviation)]

mapping_olink$Abbreviation[mapping_olink$Abbreviation == "IL-12"] <- "IL12"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "IL-18"] <- "IL18"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "IL-8"] <- "IL8"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "IL-7"] <- "IL7"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "IL-6"] <- "IL6"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "SERPINA5"] <- "SERPI5"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "SERPINA7"] <- "SERPI7"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "Protein BOC"] <- "BOC"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "THBD"] <- "TM"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "MMP-7"] <- "MMP7"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "MMP-12"] <- "MMP12"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "VEGF-A"] <- "VEGFA"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "ErbB2/HER2"] <- "ERBB2"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "ErbB3/HER3"] <- "ERBB3"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "ErbB4/HER4"] <- "ERBB4"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "FUR"] <- "FURIN"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "FasL"] <- "FASLG"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "CNTNAP2"] <- "CNTP2"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "ADAM8"] <- "ADAM-8"
mapping_olink$Abbreviation[mapping_olink$Abbreviation == "NAAA"] <- "AA"

colnames(pride_ad_prot)[!colnames(pride_ad_prot) %in% make.names(mapping_olink$Abbreviation)]

mapping_olink[!make.names(mapping_olink$Abbreviation) %in% colnames(pride_ad_prot), ]

pride_ad_unmapped <- colnames(pride_ad_prot)[
        !colnames(pride_ad_prot) %in% make.names(mapping_olink$Abbreviation)
]

sapply(
        pride_ad_unmapped,
        function(x) paste(strsplit(x, split = "\\.")[[1]][1:(length(strsplit(x, split = "\\.")[[1]]) - 1)], collapse = '.')
) %in% make.names(mapping_olink$Abbreviation)

pride_ad_nu_colnames <- mapping_olink$Uniprot[match(colnames(pride_ad_prot), make.names(mapping_olink$Abbreviation))]

pride_ad_nu_colnames[is.na(pride_ad_nu_colnames)] <- mapping_olink$Uniprot[
        match(
                sub(
                        "\\.[0-9]+$",
                        "",
                        colnames(pride_ad_prot)[is.na(pride_ad_nu_colnames)]
                ),
                make.names(mapping_olink$Abbreviation)
        )
]

nu_colnames_dups <- pride_ad_nu_colnames[duplicated(pride_ad_nu_colnames)]

pride_ad_prot_dups <- pride_ad_prot[, pride_ad_nu_colnames %in% nu_colnames_dups]

pride_ad_prot_dups_long <- t(pride_ad_prot_dups) %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(
                cols = -protein,
                names_to = "sample",
                values_to = "value"
        ) %>% as.data.frame()


ggplot(pride_ad_prot_dups_long, mapping = aes(y = protein, x = value)) +
        geom_boxplot()

pride_ad_prot_dedup <- t(pride_ad_prot)

rownames(pride_ad_prot_dedup) <- pride_ad_nu_colnames

pride_ad_prot_dedup <- t(
        sapply(
                unique(rownames(pride_ad_prot_dedup)),
                function(x) colMeans(
                        pride_ad_prot_dedup[x, , drop = FALSE], na.rm = TRUE
                )
        )
)

################################################################################
# Look at the intersection of the parsed dataframes
################################################################################

# Homogenize metadata files and match them to protein matrices, and check age
# cor distribution
################################################################################

# Emory diversity
rownames(emory_div_prot) <- sapply(
        rownames(emory_div_prot), function(x) strsplit(x, split = "\\|")[[1]][2]
)


emory_div_mdat_ord <- emory_div_mdat[
        match(colnames(emory_div_prot), emory_div_mdat$TMT.ID),
]

dim(emory_div_prot) # 1840 proteins

plot(density(cor(t(emory_div_prot), as.numeric(emory_div_mdat_ord$Age), use = "pairwise.complete.obs"), na.rm = TRUE))
mean(cor(t(emory_div_prot), as.numeric(emory_div_mdat_ord$Age), use = "pairwise.complete.obs"))
sd(cor(t(emory_div_prot), as.numeric(emory_div_mdat_ord$Age), use = "pairwise.complete.obs")) # SD ~0.15

write.csv(emory_div_prot, file = file.path(out_dir, "emory_diversity_prot_parsed.csv"))
write.csv(emory_div_mdat_ord, file = file.path(out_dir, "emory_diversity_mdat_parsed.csv"))

# Bader

dim(bader_timstof_no_dups) # 2412 proteins

plot(density(cor(t(log2(bader_timstof_no_dups)), bader_mdat$Age, use = "pairwise.complete.obs")))
mean(cor(t(log2(bader_timstof_no_dups)), bader_mdat$Age, use = "pairwise.complete.obs"))
sd(cor(t(log2(bader_timstof_no_dups)), bader_mdat$Age, use = "pairwise.complete.obs")) # SD ~0.15

write.csv(bader_timstof_no_dups, file = file.path(out_dir, "bader_prot_parsed.csv"))
write.csv(bader_mdat, file = file.path(out_dir, "bader_mdat_parsed.csv"))

# Emory 3
dim(emory_3_csf_prot) # 2229 proteins

emory_3_csf_mdat_ord <- emory_3_csf_mdat[match(colnames(emory_3_csf_prot), emory_3_csf_mdat$GUID), ]

plot(density(cor(t(emory_3_csf_prot), emory_3_csf_mdat_ord$Age, use = "pairwise.complete.obs")))
mean(cor(t(emory_3_csf_prot), emory_3_csf_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(cor(t(emory_3_csf_prot), emory_3_csf_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.15

colnames(emory_3_csf_prot_no_isos) <- make.names(colnames(emory_3_csf_prot_no_isos))
write.csv(emory_3_csf_prot_no_isos, file = file.path(out_dir, "emory_3_parsed.csv"))
write.csv(emory_3_csf_mdat_ord, file = file.path(out_dir, "emory_3_mdat_parsed.csv"))

# syn20821165_csf_disc_coh_fracTMT
dim(syn20821165_csf_disc_coh_fracTMT_prot_no_isos) # 2828 proteins

syn20821165_csf_disc_coh_fracTMT_mdat_ord <- syn20821165_csf_disc_coh_fracTMT_mdat[
        match(
                colnames(syn20821165_csf_disc_coh_fracTMT_prot_no_isos),
                sapply(
                        syn20821165_csf_disc_coh_fracTMT_mdat$SampleID,
                        function(x) strsplit(x, split = "_")[[1]][3]
                )
        ),
]

syn20821165_csf_disc_coh_fracTMT_mdat_ord$sample <- sapply(
        syn20821165_csf_disc_coh_fracTMT_mdat_ord$SampleID,
        function(x) strsplit(x, split = "_")[[1]][3]
)

syn20821165_csf_disc_coh_fracTMT_mdat_ord$Age <- as.numeric(
        gsub(
                '+',
                '',
                syn20821165_csf_disc_coh_fracTMT_mdat_ord$Age,
                fixed = TRUE
        )
)



plot(density(cor(t(syn20821165_csf_disc_coh_fracTMT_prot_no_isos), as.numeric(syn20821165_csf_disc_coh_fracTMT_mdat_ord$Age), use = "pairwise.complete.obs"), na.rm = TRUE))
mean(cor(t(syn20821165_csf_disc_coh_fracTMT_prot_no_isos), syn20821165_csf_disc_coh_fracTMT_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(
        cor(t(syn20821165_csf_disc_coh_fracTMT_prot_no_isos), as.numeric(syn20821165_csf_disc_coh_fracTMT_mdat_ord$Age), use = "pairwise.complete.obs")
) # SD ~0.05 --> Might have been regressed out?

# syn20821165_coh_1
dim(syn20821165_coh_1_prot_no_isos) # 787 proteins

syn20821165_coh_1_mdat

colnames(syn20821165_coh_1_prot_no_isos)

syn20821165_coh_1_mdat$sample <- sapply(
        syn20821165_coh_1_mdat$batch.channel,
        function(x) strsplit(x, split = '_')[[1]][3]
)

syn20821165_coh_1_mdat_ord <- syn20821165_coh_1_mdat[
        match(colnames(syn20821165_coh_1_prot_no_isos), syn20821165_coh_1_mdat$sample),
]

syn20821165_coh_1_mdat_ord$Age <- as.numeric(syn20821165_coh_1_mdat_ord$Age)

plot(density(cor(t(syn20821165_coh_1_prot_no_isos), syn20821165_coh_1_mdat_ord$Age, use = "pairwise.complete.obs")))
mean(cor(t(syn20821165_coh_1_prot_no_isos), syn20821165_coh_1_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(cor(t(syn20821165_coh_1_prot_no_isos), syn20821165_coh_1_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.08


# syn20821165_coh_2
dim(syn20821165_coh_2_prot_no_isos) # 528 proteins

syn20821165_coh_2_mdat

colnames(syn20821165_coh_2_prot_no_isos)

syn20821165_coh_2_mdat$sample <- rownames(syn20821165_coh_2_mdat)

syn20821165_coh_2_mdat_ord <- syn20821165_coh_2_mdat[
        match(colnames(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat$sample),
]

syn20821165_coh_2_mdat_ord$Age <- as.numeric(syn20821165_coh_2_mdat_ord$Age)

plot(density(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")))
mean(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.10

# syn20821165_coh_3
dim(syn20821165_coh_3_prot_no_isos) # 692 proteins

syn20821165_coh_3_mdat

colnames(syn20821165_coh_3_prot_no_isos)

syn20821165_coh_3_mdat$sample <- syn20821165_coh_3_mdat$SubjectID

syn20821165_coh_3_mdat_ord <- syn20821165_coh_3_mdat[
        match(colnames(syn20821165_coh_3_prot_no_isos), syn20821165_coh_3_mdat$sample),
]

# No age info on this metadata file.

#syn20821165_coh_3_mdat_ord$Age <- as.numeric(syn20821165_coh_3_mdat_ord$Age)


#plot(density(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")))
#mean(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs"))
#sd(cor(t(syn20821165_coh_2_prot_no_isos), syn20821165_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.10

# syn20821165_coh_4
dim(syn20821165_coh_4_prot_no_isos) # 472 proteins

syn20821165_coh_4_mdat

colnames(syn20821165_coh_4_prot_no_isos)

syn20821165_coh_4_mdat$sample <- rownames(syn20821165_coh_4_mdat)

syn20821165_coh_4_mdat_ord <- syn20821165_coh_4_mdat[
        match(colnames(syn20821165_coh_4_prot_no_isos), syn20821165_coh_4_mdat$sample),
]
syn20821165_coh_4_mdat_ord$Age <- gsub('+', '', syn20821165_coh_4_mdat_ord$Age, fixed = TRUE)
syn20821165_coh_4_mdat_ord$Age <- round(as.numeric(syn20821165_coh_4_mdat_ord$Age))

plot(density(cor(t(syn20821165_coh_4_prot_no_isos), syn20821165_coh_4_mdat_ord$Age, use = "pairwise.complete.obs")))
mean(cor(t(syn20821165_coh_4_prot_no_isos), syn20821165_coh_4_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(cor(t(syn20821165_coh_4_prot_no_isos), syn20821165_coh_4_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.10

# syn20933797_coh_2
dim(syn20933797_coh_2_prot_no_isos) # 787 proteins

syn20933797_coh_2_mdat

colnames(syn20933797_coh_2_prot_no_isos)

syn20933797_coh_2_mdat$sample <- rownames(syn20933797_coh_2_mdat)

syn20933797_coh_2_mdat_ord <- syn20933797_coh_2_mdat[
        match(colnames(syn20933797_coh_2_prot_no_isos), syn20933797_coh_2_mdat$sample),
]

plot(density(cor(t(syn20933797_coh_2_prot_no_isos), syn20933797_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")))
mean(cor(t(syn20933797_coh_2_prot_no_isos), syn20933797_coh_2_mdat_ord$Age, use = "pairwise.complete.obs"))
sd(cor(t(syn20933797_coh_2_prot_no_isos), syn20933797_coh_2_mdat_ord$Age, use = "pairwise.complete.obs")) # SD ~0.14




# Merge expression matrices

table(emory_div_mdat_ord$Diagnosis)                                # N Controls: 69
table(bader_mdat$Diagnosis_group)                              # N Controls: 599
table(emory_3_csf_mdat_ord$Diagnosis)                              # N Controls: 18
table(syn20821165_csf_disc_coh_fracTMT_mdat_ord$Diagnosis)         # N Controls: 19
table(syn20821165_coh_1_mdat_ord$Group)                            # N Controls: 63
table(syn20821165_coh_2_mdat_ord$Group)                            # N Controls: 147 # This one looks the same as syn20933797_coh_1
table(syn20821165_coh_3_mdat_ord$Group)                            # N Controls: 78
table(syn20821165_coh_4_mdat_ord$Diagnosis)                        # N Controls: 18
table(syn20933797_coh_1_mdat$Group)                            # N Controls; 147 # This one looks the same as syn20821165_coh_2
table(syn20933797_coh_2_mdat_ord$ClinicalGroup)                    # N Controls: 63

head(syn20933797_coh_1_mdat) 
head(syn20821165_coh_2_mdat)
# syn20821165_coh_2 and syn20933797_coh_1 are indeed the same


knight_somascan7K_prots <- readxl::read_xlsx("/Users/guillem.santamaria/Downloads/csf_soma7k_pQTL_aptamer_info.xlsx")

adni_cruchaga <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/cruchaga_prot_info_filt.csv")
adni_emory <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/emory_protein_info.csv")

adni_cruchaga_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/cruchaga_samp_info.csv")
adni_emory_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed/emory_samp_info.csv")

table(bader_mdat$Diagnosis_group)
table(adni_emory_mdat$DIAGNOSIS)
# emory_adni 240 CN, 59 dementia and 210 MCI
table(adni_cruchaga_mdat$DIAGNOSIS)
# Cruchaga adni 173 CN, 146 dementia and 418 MCI
table(emory_div_mdat_ord$Diagnosis)
table(emory_3_csf_mdat_ord$Diagnosis)

hist(emory_div_mdat_ord$Age)

mean(emory_div_mdat_ord$Age, na.rm = TRUE)
hist(emory_3_csf_mdat_ord$Age)
hist(bader_mdat$Age)
hist(adni_cruchaga_mdat$age)
mean(adni_cruchaga_mdat$age, na.rm = TRUE)
hist(adni_emory_mdat$age)
mean(adni_emory_mdat$age, na.rm = TRUE)


intersect_prots <- rownames(bader_timstof_no_dups) %>%
        intersect(adni_emory$uniprot_emory) %>%
        intersect(adni_cruchaga$UniProt) %>%
        intersect(mapping$Accession) %>%
        intersect(rownames(emory_3_csf_prot_no_isos)) # %>%
        #intersect(knight_somascan7K_prots$UniProt)
        #intersect(rownames(syn20821165_csf_disc_coh_fracTMT_prot_no_isos)) %>%
        #intersect(rownames(syn20821165_coh_1_prot_no_isos))# %>%
        #intersect(rownames(pride_ad_prot_dedup))
        #intersect(rownames(syn20821165_coh_2_prot_no_isos)) %>%
        #intersect(rownames(syn20821165_coh_3_prot_no_isos)) %>%
        #intersect(rownames(syn20821165_coh_4_prot_no_isos)) %>%
        #intersect(rownames(syn20933797_coh_1_prot_no_isos)) %>%
        #intersect(rownames(syn20933797_coh_2_prot_no_isos))

# Load DEG lists to see the overlap with the intersected prots
DGE_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/scRNAseq/DGE"

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


table(deg_list$cell_type)

#table(deg_list$cell_type[deg_list$gene %in% prot_info$ensembl_gene_id])


mapping_intersection <- mapping[match(intersect_prots, mapping$Accession), ]

mapping_intersection$ensembl_gene_id <- gsub(
        "\\..*",
        "",
        mapping_intersection$ensembl_gene_id
)

mapping_intersection$cell_type <- deg_list$cell_type[match(mapping_intersection$ensembl_gene_id, deg_list$gene)]

table(mapping_intersection$cell_type)


# Let's merge the following datasets, as they have a decent intersection

combined_prot_dat <- cbind.data.frame(
        emory_div_prot[intersect_prots, ],
        bader_timstof_no_dups[intersect_prots, ],
        syn20821165_csf_disc_coh_fracTMT_prot_no_isos[intersect_prots, ],
        emory_3_csf_prot_no_isos[intersect_prots, ]
)
emory_div_mdat_ord_toMerge <- emory_div_mdat_ord[,
        c("TMT.ID", "Diagnosis", "Age", "Sex", "Race", "Study")
]
bader_mdat_toMerge <- bader_mdat[, c("Run_timsTOFMeasurement", "Diagnosis_group", "Age", "Sex", "Sample_preparation_batch")]
bader_mdat_toMerge$Race <- NA
syn20821165_coh_1_mdat_ord[, c("sample", "Group", "Age", "Sex", "Race")]

colnames(emory_div_prot)

#combined_mdat <- 




# Bader dataset has been normalized with DirectLFQ
plot(density(as.matrix(bader_timstof_no_dups), na.rm = TRUE))
plot(density(as.matrix(emory_div_prot), na.rm = TRUE)) # Maybe normalized based on pool sample? 
plot(density(as.matrix(emory_3_csf_prot_no_isos), na.rm = TRUE)) # Emory 3 has some sort of normalization.
# As per the paper (https://link.springer.com/article/10.1186/s13195-022-01113-5#Sec10),
# this file (which is cleanDat.MS.CSFneat) has underwent batch correction by dividing abundances for each
# protein within each batch by the global internal standard, and then proteins with more than 75% of NAs.
plot(density(syn20821165_csf_disc_coh_fracTMT_prot_no_isos, na.rm = TRUE)) # This has some sort of 
plot(density(syn20821165_coh_1_prot_no_isos, na.rm = TRUE))
plot(density(as.matrix(syn20821165_coh_2_prot_no_isos), na.rm = TRUE))
plot(density(as.matrix(syn20821165_coh_3_prot_no_isos), na.rm = TRUE))
plot(density(as.matrix(syn20821165_coh_4_prot_no_isos), na.rm = TRUE))
plot(density(as.matrix(syn20933797_coh_1_prot_no_isos), na.rm = TRUE))
plot(density(as.matrix(syn20933797_coh_2_prot_no_isos), na.rm = TRUE))





idx <- sample(ncol(bader_timstof_no_dups), 39)
cor_sub <- cor(t(bader_timstof_no_dups[, idx]), bader_mdat$Age[idx], use = "pairwise.complete.obs")
sd(cor_sub, na.rm = TRUE)

dim(syn20821165_csf_disc_coh_fracTMT_mdat)

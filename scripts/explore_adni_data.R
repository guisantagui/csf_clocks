################################################################################
# Explore ADNI data
################################################################################
library(Hmisc)
library(dplyr)
library(stringr)
if (!require("assertr", quietly = TRUE)) {
        install.packages("assertr")
}
if (!require("bookdown", quietly = TRUE)) {
        install.packages("bookdown")
}

if (!require("ADNIMERGE2", quietly = TRUE)) {
        install.packages("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/ADNIMERGE2.tar.gz", repos = NULL, type = "source")
}

library(ADNIMERGE2)

out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results/adni_parsed"

DATADIC[DATADIC$TBLNAME == "DXSUM" & DATADIC$FLDNAME == "VISCODE", ]$MAPPING_NOTES

table(DXSUM$VISCODE)

dx_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/DXSUM_04Jun2026.csv"
data_dic_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/DATADIC_03Jun2026.csv"
demog_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/PTDEMOG_04Jun2026.csv"
prot_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/ADNI_CSFMRM_CONSOLIDATED_DATA_10Apr2014.xlsx"

nulisa_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/BSHRI_PLA_CSF_NULISA_CNS_05Jun2026.csv"

cruchaga_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/CruchagaLab_CSF_SOMAscan7k_Protein_matrix_postQC_20230620.csv"
cruchaga_protinfo_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/ADNI_Cruchaga_lab_CSF_SOMAscan7k_analyte_information_20_06_2023.csv"

emory_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/IDA/ADNI/EMORY_CSF_TMT_MS.csv"


dx <- read.csv(dx_f)
data_dic <- read.csv(data_dic_f)
demog <- read.csv(demog_f)

View(DATADIC)

# ADAE has info about age, sex and race. ADSL has this info too
dim(ADAE)

# DXSUM has info about the disease
table(DXSUM$DIAGNOSIS)


utils_file_path <- system.file(
        "analysis-dataset-utils.R",
        package = "ADNIMERGE2",
        mustWork = TRUE
)        
source(utils_file_path)


dim(readxl::read_xlsx(prot_f, sheet = 10))


cruchaga <- read.csv(cruchaga_f)
emory <- read.csv(emory_f)
nulisa <- read.csv(nulisa_f)
cruchaga_protinfo <- read.csv(cruchaga_protinfo_f)


cruchaga[1:10, 1:10]

match(cruchaga$RID, ADAE$SUBJID)

table(dx$DIAGNOSIS[match(cruchaga$RID, dx$RID)])
table(dx$DIAGNOSIS)

table(DXSUM$DIAGNOSIS[match(cruchaga$RID, DXSUM$RID)])
table(DXSUM$DIAGNOSIS)

table(DXSUM$DIAGNOSIS[match(unique(nulisa$RID), DXSUM$RID)])

unique(nulisa$UniProtID)

intersect(emory$RID, cruchaga$RID)


table(DXSUM$DIAGNOSIS[match(emory$RID, DXSUM$RID)])


cruchaga_parsed <- cruchaga
#rownames(cruchaga_parsed) <- cruchaga_parsed$RID


cruchaga_samp_info <- cruchaga[, 1:7]
# There are two duplicates, corresponding to a couple of individuals that have
# samples taken a few years apart
cruchaga_parsed[cruchaga_parsed$RID %in% cruchaga_parsed$RID[duplicated(cruchaga_parsed$RID)], 1:10]

rownames(cruchaga_parsed) <- cruchaga_parsed$GUSPECID
cruchaga_parsed <- cruchaga_parsed[, 8:ncol(cruchaga_parsed)]
plot(density(as.matrix(cruchaga_parsed), na.rm = TRUE))
sum(duplicated(cruchaga_protinfo$UniProt[match(colnames(cruchaga_parsed), cruchaga_protinfo$Analytes)]))

cruchaga_uniprot <- cruchaga_protinfo$UniProt[match(colnames(cruchaga_parsed), cruchaga_protinfo$Analytes)]
cruchaga_uniprot_dups <- cruchaga_uniprot[duplicated(cruchaga_uniprot)]
length(cruchaga_uniprot[cruchaga_uniprot %in% cruchaga_uniprot_dups])

# Let's see what do they are in protein info
cruchaga_protinfo_dups <- cruchaga_protinfo[cruchaga_protinfo$UniProt %in% cruchaga_uniprot_dups, ]
cruchaga_protinfo_dups <- cruchaga_protinfo_dups[order(cruchaga_protinfo_dups$UniProt), ]

cruchaga_protinfo_dups[cruchaga_protinfo_dups$UniProt == "P24821", ]

# ~1500 proteins in the Cruchaga dataset have duplicated UniProt IDs


which.max(table(cruchaga_protinfo_dups$UniProt))

# They are different aptamers detecting the same protein in different places.
# Let's just keep the one with better signal

# Remove probes with more than 20% NAs
cruchaga_parsed <- cruchaga_parsed[, apply(cruchaga_parsed, 2, function(x) sum(is.na(x))/length(x)) < .2]

sum(is.na(cruchaga_parsed[, cruchaga_protinfo$UniProt == cruchaga_protinfo$UniProt[1], drop = FALSE]))/nrow(cruchaga_parsed)

cruchaga_parsed_log2 <- log2(cruchaga_parsed)

cruchaga_protinfo_filt <- cruchaga_protinfo[match(colnames(cruchaga_parsed_log2), cruchaga_protinfo$Analytes), ]


keep_prots <- unlist(
        sapply(
                unique(cruchaga_uniprot),
                function(x) names(
                        which.max(
                                apply(
                                        cruchaga_parsed_log2[, cruchaga_protinfo_filt$Analytes[cruchaga_protinfo_filt$UniProt == x], drop = FALSE],
                                        2,
                                        var,
                                        na.rm = TRUE
                                )
                        )
                )
        )
)

cruchaga_parsed_log2 <- cruchaga_parsed_log2[, keep_prots]
cruchaga_parsed <- cruchaga_parsed[, keep_prots]

plot(density(as.matrix(cruchaga_parsed_log2), na.rm = TRUE))

apply(log2(cruchaga_parsed[, cruchaga_uniprot == cruchaga_uniprot_dups[1]]), 2, var, na.rm = TRUE)
apply(log2(cruchaga_parsed[, cruchaga_uniprot == cruchaga_uniprot_dups[1]]), 2, mean, na.rm = TRUE)
apply(log2(cruchaga_parsed[, cruchaga_uniprot == cruchaga_uniprot_dups[1]]), 2, function(x) sum(is.na(x))/length(x))


cruchaga_samp_info <- cbind.data.frame(
        cruchaga_samp_info,
        as.data.frame(
                DXSUM[
                        match(cruchaga_samp_info$RID, DXSUM$RID),
                        c(
                                "DIAGNOSIS",
                                "VISCODE",
                                "VISCODE2",
                                "DXMDES",
                                "DXMPTR1",
                                "DXMPTR2",
                                "DXMPTR3",
                                "DXMPTR4", "DXMPTR5", "DXMPTR6", "DXMDUE", "DXAD", "DXPARK", "DXDEP"
                        )
                ]
        )
)

cruchaga_samp_info$age <- ADSL$AGE[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$sex <- ADSL$SEX[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$race <- ADSL$RACE[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$ethnicity <- ADSL$ETHNIC[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$orig_prot <- ADSL$ORIGPROT[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$marrystat <- ADSL$APOE[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$educ <- ADSL$EDUC[match(cruchaga_samp_info$RID, ADSL$SUBJID)]
cruchaga_samp_info$apoe <- ADSL$APOE[match(cruchaga_samp_info$RID, ADSL$SUBJID)]

cruchaga_protinfo_filt <- cruchaga_protinfo_filt[match(colnames(cruchaga_parsed_log2), cruchaga_protinfo_filt$Analyte), ]

colnames(cruchaga_parsed) <- cruchaga_protinfo_filt$UniProt[match(colnames(cruchaga_parsed), cruchaga_protinfo_filt$Analytes)]

write.csv(cruchaga_parsed, file = file.path(out_dir, "cruchaga_csf_parsed.csv"))
write.csv(cruchaga_samp_info, file = file.path(out_dir, "cruchaga_samp_info.csv"))
write.csv(cruchaga_protinfo_filt, file = file.path(out_dir, "cruchaga_prot_info_filt.csv"))

# Parse Emory samples

emory_parsed <- emory[!is.na(emory$RID), ]
emory_parsed <- emory_parsed[!emory_parsed$RID %in% cruchaga_samp_info$RID, ]
emory_parsed_samp_info <- emory_parsed[, 1:7]
intersect(emory$RID, cruchaga_samp_info$RID)

emory_parsed_samp_info <- cbind.data.frame(
        emory_parsed_samp_info,
        as.data.frame(
                DXSUM[
                        match(emory_parsed_samp_info$RID, DXSUM$RID),
                        c(
                                "DIAGNOSIS",
                                "VISCODE",
                                "VISCODE2",
                                "DXMDES",
                                "DXMPTR1",
                                "DXMPTR2",
                                "DXMPTR3",
                                "DXMPTR4", "DXMPTR5", "DXMPTR6", "DXMDUE", "DXAD", "DXPARK", "DXDEP"
                        )
                ]
        )
)


emory_parsed_samp_info$age <- ADSL$AGE[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$sex <- ADSL$SEX[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$race <- ADSL$RACE[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$ethnicity <- ADSL$ETHNIC[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$orig_prot <- ADSL$ORIGPROT[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$marrystat <- ADSL$APOE[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$educ <- ADSL$EDUC[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]
emory_parsed_samp_info$apoe <- ADSL$APOE[match(emory_parsed_samp_info$RID, ADSL$SUBJID)]

rownames(emory_parsed) <- emory_parsed_samp_info$ProteomicID
emory_parsed <- emory_parsed[, 8:ncol(emory_parsed)]
emory_prot_info <- data.frame(
        symbol_emory = gsub("\\_.*", "", colnames(emory_parsed)),
        uniprot_emory = gsub(".*_", "", colnames(emory_parsed))
)

intersect(emory_prot_info$uniprot_emory, cruchaga_protinfo_filt$UniProt)

emory_prot_info$uniprot_emory[!emory_prot_info$uniprot_emory %in% cruchaga_protinfo_filt$UniProt]

library(biomaRt)

emory_mart_f <- file.path(out_dir, "emory_mapping.csv")
if (!file.exists(emory_mart_f)){
        if (!exists("mart")){
                mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
        emory_mapping <- getBM(
                attributes = c("uniprotswissprot", "uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id"),
                filters = "uniprotswissprot",
                values = emory_prot_info$uniprot_emory,
                mart = mart
        )
        write.csv(emory_mapping, file = emory_mart_f)
}else{
        emory_mapping <- read.csv(emory_mart_f)
}

colnames(emory_parsed) <- emory_prot_info$uniprot_emory

emory_prot_info <- cbind.data.frame(
        emory_prot_info,
        emory_mapping[
                match(
                        emory_prot_info$uniprot_emory,
                        emory_mapping$uniprotswissprot
                ),
                c("uniprotsptrembl", "hgnc_symbol", "ensembl_gene_id")
        ]
)

gsub("\\_.*", "", colnames(emory_parsed))[!gsub("\\_.*", "", colnames(emory_parsed)) %in% cruchaga_protinfo_filt$EntrezGeneSymbol]



write.csv(
        emory_parsed_samp_info,
        file = file.path(out_dir, "emory_samp_info.csv")
)

write.csv(
        emory_parsed,
        file = file.path(out_dir, "emory_csf_parsed.csv")
)

write.csv(
        emory_prot_info,
        file = file.path(out_dir, "emory_protein_info.csv")
)


table(emory_parsed_samp_info$DIAGNOSIS)
table(cruchaga_samp_info$DIAGNOSIS)

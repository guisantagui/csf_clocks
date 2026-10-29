library(data.table)
library(biomaRt)

prot_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/bader/DA-F08.4_-SEC-pass_v06Sc_ion_LibPGQVal1perc_precdLFQv03s01_prot_preprSc04.tsv"
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/results"

read_table_fast <- function(f, row.names = NULL, sep = "auto", sep2 = "auto"){
        df <- data.frame(fread(f, sep = sep, sep2 = sep2))
        if (!is.null(row.names)){
                rownames(df) <- df[, row.names]
                df <- df[, -row.names]
        }
        return(df)       
}
prot <- read_table_fast(prot_file, sep = "\t")

human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
        attributes = c("hgnc_symbol", "ensembl_gene_id", "description"),
        filters = "hgnc_symbol",
        values = prot$Genes,
        mart = human
)

write.table(mapping, file = file.path(out_dir, "id_mapping.tsv"), row.names = FALSE, col.names = TRUE)

prot[1:10, 1:10]

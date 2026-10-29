synID_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/CSF Studies - SynIDs.csv"
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/"

synIDs <- read.csv(synID_f)$Synapse.ID

synID_df <- data.frame(dataset = "CSF", synID = synIDs)

write.table(synID_df, file = file.path(out_dir, "synIDs_parsed.tsv"), sep = '\t', quote = F, col.names = F, row.names = F)

ori_file <- file.path(out_dir, "synIDs_parsed.tsv")

setwd("/Users/guillem.santamaria/Documents/postdoc/bulkRNAseq_pipe/scripts")

visited_synIDs <- c()

getChildren_rec <- function(syn_file, out_dir){
        file_dir <- dirname(syn_file)
        #synIDs <- read.table(ori_file)$V2
        if (file.exists(syn_file) && file.info(syn_file)$size != 0) {
                synIDs <- read.table(syn_file, sep = '\t')$V2
                for(synID in synIDs){
                        system(sprintf("Rscript get_synID_children.R --synID %s --out_dir %s", synID, out_dir))
                        getChildren_rec(
                                file.path(file_dir, sprintf("%s_children.tsv", synID)),
                                out_dir
                        )
                }
        }else{
                endpoint <- gsub("_children.tsv", "", basename(syn_file))
                print(sprintf("Endpoint: %s", endpoint))
        }
}
read.table("/Users/guillem.santamaria/Documents/postdoc/comput/csf_clocks/data/csf_proteomics/synapse/syn44132374_children.tsv", sep = '\t')
getChildren_rec(ori_file, out_dir = out_dir)

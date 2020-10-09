source("Scripts/plink_tools.R")
library(dplyr)

find_missing <- function(){
    total_df <- data.frame()
    for (filename in Sys.glob("*.frqx")) {
        current_df <- read_delim(filename, "\t")
        
        ## Load in the corresponding Hardy Weinberg data and add the P-value column to the table
        hwe_df <- read_table2(Sys.glob(paste0("*", get_root_name(filename), ".hwe")))
        current_df <- cbind(current_df, "HWE P" = hwe_df$P)
        
        ## Subset the data to only keep variants with two or less homozygous cases
        sub_df <- current_df[current_df$`C(HOM A1)` <= 2 | current_df$`C(HOM A2)` <= 2,]
        if (nrow(sub_df) == 0){
            next()
        }
        
        ## Lookup the remaining variants in the VCF file and append the INFO column to the table
        vcf_df <- read_vcf(Sys.glob(paste0("*", get_root_name(filename), ".*vcf")))
        vcf_df <- merge(sub_df, vcf_df, by.x = "SNP", by.y = "POS")
        sub_df <- cbind(sub_df, "INFO" = vcf_df$INFO)
        
        ## Add the gene name to the dataframe and append to the main file
        Gene <- c(get_root_name(filename, TRUE))
        sub_df <- cbind(Gene, sub_df)
        total_df <- rbind(total_df, sub_df)
    }
    total_df
    write.csv(total_df, "top_variants.csv", row.names = FALSE)
}

main <- function(){
    wd <- getwd()
    setwd("Data/")
    
    #make_tfiles()
    #fix_tped()
    #run_cmd("--freqx --hardy")
    find_missing()
    
    setwd(wd)
}

summarize <- function() {
    top_variants <- read_csv("Data/top_variants.csv")
    
    Counts <- c("Total", "HWE <0.05", "HWE <0.01", "Modifier", "Low", "Moderate", "High")
    total <- table(top_variants$Gene)
    HWE5 <- table(top_variants[top_variants$`HWE P`<=0.05,]$Gene)
    HWE1 <- table(top_variants[top_variants$`HWE P`<=0.01,]$Gene)
    mod <- table(top_variants[grepl("MODIFIER", top_variants$INFO, fixed = TRUE),]$Gene) 
    low <- table(top_variants[grepl("LOW", top_variants$INFO, fixed = TRUE),]$Gene) 
    mid <- table(top_variants[grepl("MODERATE", top_variants$INFO, fixed = TRUE),]$Gene) 
    high <- table(top_variants[grepl("HIGH", top_variants$INFO, fixed = TRUE),]$Gene) 
    
    comb <- bind_rows(total, HWE5, HWE1, mod, low, mid, high)
    comb <- cbind(counts, comb)
    
    write.csv(comb, "summary.csv", row.names = FALSE, na = "0")
}

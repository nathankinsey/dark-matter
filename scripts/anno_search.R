library(readr)
library(dplyr)
library(tidyr)
library(stringr)

init <- function() {
    # Download files for looking up CpG islands, regulatory regions, and exon splice sites
    dir.create("Data/")
    download.file("http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/cpgIslandExt.txt.gz",
                  "Data/cpgIsland.txt.gz")
    download.file("http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/oreganno.txt.gz",
                  "Data/oreganno.txt.gz")
    download.file("http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/oregannoAttr.txt.gz",
                  "Data/oregannoAttr.txt.gz")
    download.file('http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "clfamiliaris_gene_ensembl" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "exon_chrom_start" /><Attribute name = "exon_chrom_end" /><Attribute name = "ensembl_exon_id" /></Dataset></Query>',
                  "Data/ensemblExons.txt")
    
    # Combine all data sources into a master annotation reference file
    cpgIsland <- read_delim("Data/cpgIsland.txt.gz", "\t", escape_double = FALSE,
                            col_names = c("bin", "chrom", "bpStart", "bpEnd", "name",
                                          "length", "cpgNum", "gcNum", "perCpg", "perGc",
                                          "obsExp"), 
                            trim_ws = TRUE)
    oreganno <- read_delim("Data/oreganno.txt.gz", "\t", escape_double = FALSE,
                           col_names = c("bin", "chrom", "bpStart", "bpEnd", "id",
                                         "strand", "name"), 
                           trim_ws = TRUE)
    oregannoAttr <- read_delim("Data/oregannoAttr.txt.gz", "\t", escape_double = FALSE,
                           col_names = c("id", "attribute", "attrVal"), trim_ws = TRUE)
    ensemblExons <- read_delim("Data/ensemblExons.txt", "\t", escape_double = FALSE,
                               col_names = c("chrom", "bpStart", "bpEnd", "exonName"), 
                               trim_ws = TRUE)
    
    cpgIsland <- cpgIsland %>%
        mutate(notes = "CpG Island") %>%
        select(chrom:bpEnd, notes) %>%
        mutate(chrom = str_remove(chrom, "chr"))
    
    oregannoAttr <- spread(oregannoAttr, attribute, attrVal)
    oregMaster <- inner_join(oreganno, oregannoAttr) %>%
        mutate(notes = paste0(type, " (", Gene, ")")) %>%
        mutate(chrom = str_remove(chrom, "chr")) %>%
        select(chrom:bpEnd, notes)
    rm(oreganno, oregannoAttr)
    
    ensemblExons <- ensemblExons %>%
        mutate(notes = paste0("Exon splice site (", exonName, ")")) %>%
        select(chrom:bpEnd, notes)
    
    annoRef <- bind_rows(cpgIsland, oregMaster, ensemblExons) %>%
        distinct()
    
    write_csv(annoRef, "annoRef.csv", col_names = TRUE)
}

# Note that chrom will only match the X chromosome if given a capital "X" as input.
# Entering the corresponding number will not yield results, as that is how it is encoded in UCSC and Ensembl.
annoSearch <- function(chr, bp, id = NULL, searchWidth = 500, export = NULL) {
    annoRef <- read_csv("annoRef.csv", col_types = cols(chrom = col_character()))
    
    return_DF <- tibble()
    for (i in 1:length(chr)) {
        sub_DF <- filter(annoRef, chr[i] == chrom,
            (bp[i] >= bpStart - searchWidth & bp[i] <= bpStart + searchWidth) |
            (bp[i] >= bpEnd - searchWidth & bp[i] <= bpEnd + searchWidth))
        
        if (!is.null(id)) {
            sub_DF <- mutate(sub_DF, closeTo = paste0(id[i], " (", chr[i], ":", bp[i], ")"))
        }
        
        return_DF <- bind_rows(return_DF, sub_DF)
    }
    
    if (is.null(export)) {
        return(return_DF)
    } else {
        write_csv(return_DF, export, col_names = TRUE)
    }
}
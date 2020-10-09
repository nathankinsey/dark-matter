library(readr)

get_root_name <- function(filename, geneName = FALSE){
    iteratable <- strsplit(filename, "")
    
    ## Find the location of the period in the filename
    i <- 1
    while(iteratable[[1]][i] != "."){
        i <- i + 1
    }
    
    ## Output only the gene name or file name, depending on the options
    if (geneName){
        substr(filename, 8, i - 1)
    } else {
        substr(filename, 1, i - 1)
    }
}

make_tfiles <- function(){
    for (filename in Sys.glob("*.vcf")) {
        system(paste0("plink --vcf ", filename, " --recode transpose tab --geno 0.1 --maf 0.05 --dog --const-fid --out ", get_root_name(filename)))
    }
}

fix_tped <- function(){
    for (filename in Sys.glob("*.tped")) {
        current_df <- read_delim(filename, "\t", col_names = FALSE)
        
        ## Copy the bp coordinate column into the variant identifier column
        current_df[,2] <- current_df[,4]
        write.table(current_df, filename, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
    }
}

run_cmd <- function(cmd, suffix = ""){
    for (filename in Sys.glob("*.tped")) {
        root_name <- get_root_name(filename)
        system(paste0("plink --tfile ", root_name, " --dog ", cmd, " --out ", root_name, suffix))
    }
}

read_vcf <- function(filename){
    ## Read the VCF line-by-line until the actual data begins and record that line # (variable by file) 
    skipNum <- 0
    con <- file(filename, "r")
    
    repeat {
        line = readLines(con, n = 1)
        
        ## Throw an error if the EOF is reached without finding the data header
        if (length(line) == 0) {
            stop("There was a problem reading the VCF ", filename)
        }
        
        if (grepl("#CHROM", line, fixed = TRUE)) {
            break()
        }
        skipNum <- skipNum + 1
    }
    close(con)
    
    read_table2(filename, skip = skipNum)
}
library(tidyverse)

dbvdcVars <- read_csv("data/dbvdcVars.csv", col_types = cols(`#CHROM` = col_character())) %>% 
    select(Chr = `#CHROM`, POS) %>%
    distinct() %>%
    drop_na()

candidateGenes <- read_csv("data/810data/ensemblCandidateData.csv") %>% 
    filter(is.na(Gene)) %>% 
    select(CHR, GeneID, bpStart, bpEnd) %>%
    distinct()

lookupGene <- function(chr, pos) {
    working <- filter(candidateGenes, chr == CHR, pos >= bpStart, pos <= bpEnd)
    if (nrow(working) == 0) {
        return(NA)
    } else {
        return(working$GeneID[1])
    }
}

uncharacterizedCandidates <- rowwise(dbvdcVars) %>% 
    transmute(gene = lookupGene(Chr, POS)) %>% 
    drop_na() %>%
    distinct()

# ----- Find Overlap Between Flagged Variants -----
library(tidyverse)

supplementalTable1 <- read_csv("~/supplementalTable1.csv") %>% 
    mutate(pos = paste0(Chromosome, ":", `Position (bp)`))

dbvdcVars <- read_csv("data/dbvdcVars.csv", col_types = cols(`#CHROM` = col_character())) %>% 
    mutate(pos = paste0(`#CHROM`, ":", POS))

supplementalTable1 <- mutate(supplementalTable1, "inDBVDC" = pos %in% dbvdcVars$pos)
supplementalTable1 <- mutate(supplementalTable1, "geneInDBVDC" = `Gene Name` %in% dbvdcVars$GENE & !is.na(`Gene Name`))

write_csv(supplementalTable1, "~/supplementalTable1.csv", na = "")

# ----- Filter out characterized genes -----
library(tidyverse)

uncharacterized <- ensemblCandidateData %>%
    filter(!(
        GeneID %in% c(
            "ENSCAFG00000000604",
            "ENSCAFG00000017072",
            "ENSCAFG00000007311",
            "ENSCAFG00000048508",
            "ENSCAFG00000018129"
        )
    )) %>% 
    filter(is.na(Gene))
library(tidyverse)

gwasData <- read_table2("data/DM_GENO1_MAF05.hwe")
arrayManifest <- read_csv("data/arrayManifest.csv")

# PPP2CB is 16:33532958-33566640

snpsInterest <- filter(arrayManifest, Chr == 16,
                       MapInfo >= 33531958, 
                       MapInfo <= 33567640)$Name

gwasFiltered <- filter(gwasData, SNP %in% snpsInterest)

# AMO

AMO687 <- read_delim("data/AMO/AMO687.tfam", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     col_types = cols_only(X2 = col_guess(), 
                                           X5 = col_guess()), trim_ws = TRUE) %>% 
    rename(dog = X2, sex = X5)

phosSNP687 <- read_table2("data/AMO/phosSNP687.vcf", 
                          col_names = FALSE, col_types = cols(X1 = col_skip(), 
                                                              X2 = col_skip(), X3 = col_skip(), 
                                                              X4 = col_skip(), X5 = col_skip(), 
                                                              X6 = col_skip(), X7 = col_skip(), 
                                                              X8 = col_skip(), X9 = col_skip()), 
                          skip = 6) %>% 
    t()
phosSNP687 <- tibble(phosSNP687[,1], phosSNP687[,2]) %>% 
    rename(ID = `phosSNP687[, 1]`, Call = `phosSNP687[, 2]`) %>% 
    separate(ID, into = c("ignore", "dog"), sep = "_") %>% 
    inner_join(AMO687)

sample(filter(phosSNP687, Call == "0/1", sex == "2")[[1]], 5)

males <- filter(AMO687, sex == 1)$dog
females <- filter(AMO687, sex == 2)$dog

# ---- Sex Differences on WGS Data ----

# Get phenotype information
sampleInfo <- fread("data/sample_info.csv") %>% 
    mutate(isWolf = BreedGroup == "Wolf") %>% 
    #select(-BreedGroup)


ppp2cb_lethal <- read_csv("ppp2cb_lethal.csv") %>% 
    separate(SampleId, into = c("ignore", "SampleId"), sep = "_")

pheno <- inner_join(sampleInfo, ppp2cb_lethal)

maleData <- filter(pheno, Sex == "m")
femaleData <- filter(pheno, Sex == "f")

table(maleData$`16:33555148`)
table(femaleData$`16:33555148`)

table(maleData$`16:33555160`)
table(femaleData$`16:33555160`)

table(maleData$`16:33555181`)
table(femaleData$`16:33555181`)

table(maleData$`16:33555184`)
table(femaleData$`16:33555184`)

table(maleData$`16:33560533`)
table(femaleData$`16:33560533`)

table(maleData$`16:33560537`)
table(femaleData$`16:33560537`)

table(maleData$`16:33560558`)
table(femaleData$`16:33560558`)

table(maleData$`16:33566330`)
table(femaleData$`16:33566330`)

table(maleData$`16:33566351`)
table(femaleData$`16:33566351`)

# ------------------------------- Use gaston to make an LD plot ------------------------------------
library(gaston)

vcf <- read.vcf("data/regions/ppp2cb_lethal_plus.vcf")
nSnps <- ncol(vcf)
ldMatrix <- LD(vcf, lim = c(1, nSnps))
ldMatrix[is.nan(ldMatrix)] <- 0

# Save LD plot to the working directory in a high quality tiff
tiff("test.tiff", units = "in", width = 5, height = 5, res = 300)
LD.plot(ldMatrix[1:nSnps,1:nSnps], snp.positions = vcf@snps$pos[1:nSnps])
dev.off()

# ---- Breed differences on WGS Data ----
sampleInfo <- fread("data/wgs/sampleInfo.csv") %>% 
    mutate(across(c("Breed", "Breed Group"), tolower))

ppp2cb_lethal <- read_csv("ppp2cb_lethal.csv") %>% 
    separate(SampleId, into = c("ignore", "Sample ID"), sep = "_") %>% 
    inner_join(sampleInfo)

total <- ppp2cb_lethal %>% 
    group_by(Breed, `16:33566330`) %>%
    summarise(n = n()) %>% 
    filter(n >= 5)

# -------------------------------------- AMO 42 WGS Sequences --------------------------------------
ppp2cb <- read_vcf("data/AMO/ppp2cb42amo.vcf")

table(ppp2cb$`16:33555148`)
table(ppp2cb$`16:33555160`)
table(ppp2cb$`16:33555181`)

hets <- filter(ppp2cb, `16:33555181` == "0/1")



# ------------------------------------ Who Female Hets in DBVDC ------------------------------------
ppp2cbC <- read_vcf("data/regions/ppp2cb_lethal_plus.vcf") %>% 
    separate(SampleId, into = c("ignore", "SampleId"), sep = "_") %>% 
    select(-ignore)

sampleInfo <- fread("data/full_pheno.csv") %>% 
    inner_join(ppp2cbC, by = c("Sample ID" = "SampleId")) %>% 
    filter(`16:33555181` == "0/1", Sex == "f")

sampleInfo <- fread("data/813DBVDCmeta.csv") %>% 
    inner_join(ppp2cbC) %>% 
    filter(str_detect(Breed, "Wolf"))

# ------------------------------ Check Sex Differences in WGS SNPs ---------------------------------
read_vcf <- function(filename) {
    # Read in a VCF and remove the unimportant data
    vcf <- fread(filename, skip = "#CHROM") %>% 
        select(!(ID:FORMAT)) %>% 
        unite("#CHROM", "POS", col = "Location", sep = ":")
    
    # Remove all information under each sample column except for called allele
    cleanAlleles <- vcf %>% 
        select(!Location) %>% 
        apply(2, str_extract, pattern = "(^[\\d.]/[\\d.])") %>% 
        as_tibble()
    vcf <- bind_cols(select(vcf, Location), cleanAlleles)
    
    # Transpose the table to tidy the data
    vcf = t(vcf)
    colnames(vcf) <- vcf[1,]
    vcf <- vcf[-1,] %>% 
        as_tibble(rownames = "SampleId")
}

sampleInfo <- fread("data/full_pheno.csv")
overlapCandidates <- read_csv("overlapCandidates.csv")

voi <- filter(overlapCandidates, Gene == "ENSCAFG00000006025" | Gene == "ENSCAFG00000049994" | 
                  Gene == "ENSCAFG00000046019") %>% 
    transmute(pos = paste0(CHR, ":", BP))
trim21 <- read_vcf("data/regions/TRIM21etc.recode.vcf") %>% 
    inner_join(sampleInfo, by = c("SampleId" = "Sample ID")) %>% 
    select(SampleId, voi$pos, Sex)
table(trim21[[10]], trim21$Sex)

for (i in 2:41) {
    print(colnames(trim21)[i])
    print(table(trim21[[i]], trim21$Sex))
}

voi <- filter(overlapCandidates, Gene == "ENSCAFG00000003646") %>% 
    transmute(pos = paste0(CHR, ":", BP))
sfrp4 <- read_vcf("data/regions/SFRP4.recode.vcf") %>% 
    inner_join(sampleInfo, by = c("SampleId" = "Sample ID")) %>% 
    select(SampleId, voi$pos, Sex)
table(sfrp4[[7]], sfrp4$Sex)


voi <- filter(overlapCandidates, Gene == "ENSCAFG00000048018") %>% 
    transmute(pos = paste0(CHR, ":", BP))
lnc <- read_vcf("data/regions/ENSCAFG00000048018.recode.vcf") %>% 
    inner_join(sampleInfo, by = c("SampleId" = "Sample ID")) %>% 
    select(SampleId, voi$pos, Sex)
table(lnc[[3]], lnc$Sex)

voi <- filter(overlapCandidates, Gene == "ENSCAFG00000030470" | Gene == "ENSCAFG00000041929" | 
                  Gene == "ENSCAFG00000049488") %>% 
    transmute(pos = paste0(CHR, ":", BP))
olf <- read_vcf("data/regions/OR52H2etc.recode.vcf") %>% 
    inner_join(sampleInfo, by = c("SampleId" = "Sample ID")) %>% 
    select(SampleId, voi$pos, Sex)

for (i in 2:8) {
    print(colnames(olf)[i])
    print(table(olf[[i]], olf$Sex))
}

# ----------------------------------- Check reads of AMO VCF ---------------------------------------
library(data.table)
library(tidyverse)
amoDat <- fread("data/AMO/BCBDBT42_chr21.eff.vcf", skip = "#CHROM")
orDat <- filter(amoDat, `#CHROM` == 21 & POS > 27089150 & POS < 27091078)
trimDat <- filter(amoDat, `#CHROM` == 21 & POS > 26775671 & POS < 26783726)
olfDat <- filter(amoDat, `#CHROM` == 21 & POS > 28539571 & POS < 28540520)
write_csv(olfDat, "OR52H2AMO.csv")

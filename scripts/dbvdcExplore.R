library(tidyverse)

dbvdcVars <- read_csv("dbvdcVars.csv", col_types = cols(`#CHROM` = col_character()))

allVarsCount <- dbvdcVars %>% 
  select(CHROM = "#CHROM", "POS") %>% 
  distinct() %>% 
  group_by(CHROM) %>% 
  transmute(allCount = n()) %>% 
  distinct() %>% 
  arrange(desc(allCount))

elCount <- dbvdcVars %>% 
  filter(IMPACT == "EL") %>% 
  select(CHROM = "#CHROM", "POS") %>% 
  distinct() %>% 
  group_by(CHROM) %>% 
  transmute(elCount = n()) %>% 
  distinct()

ddCount <- dbvdcVars %>% 
  filter(IMPACT == "DD") %>% 
  select(CHROM = "#CHROM", "POS") %>% 
  distinct() %>% 
  group_by(CHROM) %>% 
  transmute(ddCount = n()) %>% 
  distinct()

geneCount <- dbvdcVars %>% 
  select(CHROM = "#CHROM", "GENE") %>% 
  distinct() %>% 
  group_by(CHROM) %>% 
  transmute(geneCount = n()) %>% 
  distinct()

counts <- inner_join(allVarsCount, ddCount) %>% 
  inner_join(elCount) %>% 
  inner_join(geneCount)
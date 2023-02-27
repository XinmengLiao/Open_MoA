## 1)) Data Preparation

library(tidyverse)
library(dbparser)
library(Biostrings)

## String
## Screen String data with top 10% of combined scores
setwd("Datapreparation/String")
string <- read.csv('StringRaw_9606.protein.physical.links.full.v11.5.txt',sep = '',header = T)
string_exp <- string[which(string$experiments != 0 | string$experiments_transferred !=0),]
string_exp <- string_exp %>% select(protein1, protein2, combined_score)
info <- read.delim('StringRaw_9606.protein.info.v11.5.txt',sep = '\t',header = T)
info <- info[,c(1,2)]
info2 <- info
colnames(info)[1] <- 'protein1'
colnames(info2)[1] <- 'protein2'
name1 <- merge(string_exp, info, by = 'protein1')
colnames(name1)[4] <- 'name1'
name2 <- merge(name1, info2, by = 'protein2')
name2$name_1 <- 'Protein'
name2$name_2 <- 'Protein'
name2$Relationship <- 'PPI'
name2$common_score <- 1
string <- name2 %>% select(name1,name_1,preferred_name,name_2, Relationship, common_score, combined_score)
colnames(string)[c(1,3)] <- c('Vertex.A','Vertex.B')
string10 <- string[which(string$combined_score > quantile(string$combined_score, 0.9)),]
write.csv(string10,file = 'StringTop10.csv',row.names = F)


### RegNetwork 
setwd("Datapreparation")
reg <- read.table('RegnetworkRaw_human.source')
index1 <- which(grepl('^hsa-',x = reg$V1))
index2 <- which(grepl('^MIR',x = reg$V3))
index3 <- which(grepl('^hsa-',x = reg$V3))
reg_tf <- reg[-c(index1,index2,index3),]
reg_tf <- reg_tf[,c(1,3)]
colnames(reg_tf) <- c('Vertex.A','Vertex.B')
reg_tf <- reg_tf %>% 
  mutate('TF', 'Gene',1,0)
colnames(reg_tf)[3:6] <- c('name_1','name_2','common_score','combined_score')
reg_tf$Relationship <- 'Regulatory'
reg_tf <- reg_tf[,c(1,3,2,4,7,5,6)]
write.csv(reg_tf,file = '~/Desktop/RegnetworkTF.csv',row.names = F)


### DrugBank
setwd("~/Desktop/Sweden/Scorenetwork 2.0/Datapreparation/DrugBank")
read_drugbank_xml_db('full_database.xml')
cettfile <- cett(
  save_csv = T,
  save_table = F,
  csv_path = '.',
  override_csv = F,
  database_connection = NULL)
target <- read.csv('targets_polypeptides.csv')
h_target <- target[which(target$organism == 'Humans'),]
target_brief <- h_target %>% 
  select(id,name,gene_name, parent_id)
colnames(target_brief) <- c('UniProt.ID','Protein_name','Gene_name','parent_id') 
link <- read.csv('drug_target_link.csv')
link <- link %>% 
  select(DrugBank.ID, Name, UniProt.ID)
all_drugbank <- merge(link, target_brief, by = 'UniProt.ID')
all_drugbank <- all_drugbank %>% select(Name, Gene_name)
drugbank <- distinct(all_drugbank)
colnames(drugbank) <- c('Vertex.A','Vertex.B')
drugbank <- drugbank %>% 
  mutate('Drug','Gene', 'Drug_Gene',1,0)
colnames(drugbank)[3:7] <- c('name_1','name_2','Relationship','common_score','combined_score')
drugbank <- drugbank[,c(1,3,2,4,5,6,7)]
write.csv(drugbank, file = 'DrugBank.csv',row.names = F)


## 2)) Construction of Background Network
drugbank <- read.csv('DrugBank.csv')
regnetwork <- read.csv('RegnetworkTF.csv')
string <- read.csv('StringTop10.csv')
scorenetwork <- rbind(drugbank, regnetwork, string)
scorenetwork <- unique(scorenetwork)
scorenetwork <- distinct(scorenetwork)
index1 <- which(scorenetwork$Vertex.B == '')
scorenetwork <- scorenetwork[-index1,]
write.csv(scorenetwork, file = 'General_BG_Network.csv',row.names = F)


## 3)) Specify to the HepG2 Network
## here we choose the HepG2 RNA data from CCLE and HPA databases 
setwd("Datapreparation")
general_network <- read.csv('General_BG_Network.csv')
hpa_ccle <- read.csv('Cell_Line_screen/tpm_hepg2_hpa_ccle.csv',header = T)

hpa_ccle$max <-NA
for (i in 1:nrow(hpa_ccle)) {
  print(i)
  if(hpa_ccle$hpa_tpm[i] > hpa_ccle$ccle_tpm[i]){
    hpa_ccle$max[i] <- hpa_ccle$hpa_tpm[i]
  } else {
    hpa_ccle$max[i] <- hpa_ccle$ccle_tpm[i]
  }
}

hpa_ccle_max <- hpa_ccle
hpa_ccle_max <- hpa_ccle_max[which(hpa_ccle_max$max > 1),]
hpa_ccle_max <- rbind(hpa_ccle_max,c('WNT1',0,0,0,0))
hepg2_network <- 
  general_network[which(general_network$Vertex.A %in% hpa_ccle_max$symbol & general_network$Vertex.B %in% hpa_ccle_max$symbol),]
hepg2_network <- hepg2_network[which(hepg2_network$Relationship  != 'Drug_Gene'),]
drugbank <- general_network[which(general_network$Relationship == 'Drug_Gene'),]
scorenetwork <- rbind(hepg2_network, drugbank)
write.csv(scorenetwork,file = 'Hepg2_Network_top10_HpaCcle_TPMmax>1.csv',row.names = F)
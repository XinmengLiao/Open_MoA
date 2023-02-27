if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
install.packages('ggh4x')
library(ggplot2)
library(stringr)
library(tidyverse)
library(enrichR)
library(ggh4x)
library(rWikiPathways)
library(forcats)
library(RCy3)
library(dplyr)
library(ggpubr)
library(xlsx)
cytoscapePing()

### KEGG bubble plot
KEGG_raw <- read.table('WNT_KEGG_gene99.txt',sep = '\t',header = T)
KEGG_raw$Term <- str_split_fixed(KEGG_raw$Term,'\\:',2)
KEGG_raw$Term <- KEGG_raw$Term[,-1] 
KEGG <- KEGG_raw[KEGG_raw$FDR < 0.05,]
KEGG <- KEGG[grepl(x = KEGG$Term,pattern = 'signaling pathway'),]
KEGG <- KEGG[1:5,]
KEGG$log10 <- -log10(KEGG$FDR)
KEGG$Term <- fct_reorder(KEGG$Term,KEGG$PValue,.desc = T)

plot1 <- ggplot(data = KEGG, aes(x =log10,y = Term,color = FDR))+
  theme(plot.background = element_blank())+
  theme_bw(base_line_size = 1,base_rect_size = 2)+
  geom_point(stat = 'Identity',size = 6,alpha = 0.8)+
  theme(aspect.ratio = 1)+
  scale_color_gradient(low = 'red',high = 'blue')+
  ylab(label = 'KEGG Pathways')+
  xlab(label = expression(paste("-Log"[10], "(FDR)")))+
  theme(axis.title.x = element_text(size = 15,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.x = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.title.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'))+
  theme(legend.key.size = unit(0.6,'cm'),
    #legend.key.height = unit(0.3,'cm'),
        #legend.key.width = unit(1,'cm'),
        legend.text = element_text(family = 'Times',size = 10,face = 'bold'),
        legend.title = element_text(family = 'Times',size = 10,face = 'bold',vjust = 0.95),
        legend.position = 'right', legend.box.spacing = unit(0.5,'mm'))+
  scale_x_continuous(breaks = seq(5,12,1))
plot1
ggsave(plot1, filename = 'Wnt1-2.png',width = 25, height = 16,units = 'cm',dpi = 300)


#### Wikipathway analysis
pathways <- findPathwaysByText('"Wnt"') 
human.pathways <- pathways %>% 
  dplyr::filter(species == "Homo sapiens") 
url <- getPathwayInfo(c('WP363'))$url ## check whether the pathway ir right
browseURL(url)
cytoscapePing()
for (i in 3:nrow(human.tgfb1.pathways)){
  RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',human.tgfb1.pathways[i,2])) 
}
RCy3::commandsRun('wikipathways import-as-pathway id=WP560') 
toggleGraphicsDetails()

esnl <- read.table('~/Desktop/Sweden/Scorenetwork 2.0/ENSG_num.txt',header = T)
genelist <- read.table('TGFb1_gene99_TPMmax>1.txt',header = T)
colnames(genelist) <- 'symbol'
genelist2 <- left_join(genelist, esnl,by = 'symbol')
genelist2$color <- 1 
loadTableData(genelist2, data.key.column = "gene_id", table.key.column ="Ensembl")
setNodeColorBypass(node.names = c('E2F5',"SMAD2",'ATF2'),new.colors = '#FD39B8') 



### OMIM database from Enrichr 
## Enrichr 
setwd("~/Desktop/Sweden/Scorenetwork 2.0/Metformin/Metformin_2.0")
listEnrichrSites()
setEnrichrSite("Enrichr") # Enrichr 就是 Human genes，还有其他物种的可以自己根据list设置
websiteLive <- TRUE
dbs <- listEnrichrDbs()
dbs <- c("OMIM_Disease") ## choose the background database
gene <- read.csv('Metformin_gene99_TPMmax>1.txt')  ## import the gene list
if (websiteLive) {
  enriched <- enrichr(gene$OFFICAL_GENE_SYMBOL, dbs)
}

if (websiteLive) enriched[["OMIM_Disease"]]
if (websiteLive) {
  #x <- head(enriched[["OMIM_Disease"]])
  x <- enriched[["OMIM_Disease"]]
  x[,1] <- gsub("GO:", "GO_", x[,1])
  table(x)
}
write.xlsx(x,'metformin_OMIM all.xlsx')

x$log10 <- -log10(x$P.value)
colnames(x)[3] <- 'P value'
x$Term[1:nrow(x)] <- c('Type 2 diabetes mellitus',
                       'Colorectal cancer',
                       'Pancreatic cancer',
                       'Diabetes',
                       'Leukemia',
                       'Ovarian cancer')

x$Term <- fct_reorder(x$Term,x$`P value`,.desc = T)
Enrich <- ggplot(data = x, aes(x =log10,y = Term,color =`P value`))+
  theme_bw(base_size = 10, base_line_size =1,base_rect_size = 2)+
  geom_point(stat = 'Identity',size = 5,alpha = 0.8)+
  scale_color_gradient(low = 'red',high = 'blue')+
  theme(legend.text = element_text(family = 'Times',size = 11,face = 'bold'),
        legend.title = element_text(family = 'Times',size = 15,face = 'bold'),
        aspect.ratio = 1)+
  ylab(label = 'OMIM Diseases')+
  xlab(label = expression(paste("-Log"[10], "(P value)")))+
  theme(axis.title.x = element_text(size = 15,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.x = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.title.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'))+
  theme(legend.key.size = unit(0.6,'cm'),
        legend.text = element_text(family = 'Times',size = 8,face = 'bold'),
        legend.title = element_text(family = 'Times',size = 10,face = 'bold',vjust = 0.95),
        legend.position = 'right', legend.box.spacing = unit(0.5,'mm'))+
  scale_x_continuous(breaks = seq(1,2,0.1))
Enrich
met_OMIM <- Enrich

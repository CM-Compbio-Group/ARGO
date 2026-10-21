#-------------------------------------
## run enrichment analysis and plot results
#-------------------------------------
library(fgsea)
library(msigdbr)

## Hallmark genesets
h_gene_sets = msigdbr(db_species = 'MM', species = "Mus musculus", collection = "MH")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
names(msigdbr_list) <- str_sub(names(msigdbr_list), start = 10)
# msigdbr_list$E2F_TARGETS

## GO genesets
go_gene_sets = msigdbr(db_species = 'MM', species = "Mus musculus", collection = "M5", subcollection = 'GO:BP')
go_list = split(x = go_gene_sets$gene_symbol, f = go_gene_sets$gs_name)
names(go_list) <- str_sub(names(go_list), start = 6)

## KEGG genesets
reactome_gene_sets = msigdbr(db_species = 'MM', species = "Mus musculus", collection = "M2", subcollection = 'CP:REACTOME')
reactome_list = split(x = reactome_gene_sets$gene_symbol, f = reactome_gene_sets$gs_name)
names(reactome_list) <- str_sub(names(reactome_list), start = 10)

#------------------------------------------
## get list of diff genes
#------------------------------------------
ranks <- funGP_res$stats 
names(ranks) <- funGP_res$gene
ranks <- ranks[!is.na(ranks)]
plot(sort(ranks))

library(extrafont)
hallmark_gs <- fgsea(pathways = msigdbr_list, stats = ranks, minSize = 10, maxSize = 500)
make_enrichment_dotplot(hallmark_gs, 'Hallmark', wd=7, ht=7)

go_gs <- fgsea(pathways = go_list, stats = ranks, minSize = 10, maxSize = 500)
make_enrichment_dotplot(go_gs, 'GO-BP', 
                       wd = 7,ht = 7)

reactome_gs <- fgsea(pathways = reactome_list, stats = ranks, minSize = 10, maxSize = 500)
make_enrichment_dotplot(reactome_gs, 'REACTOME', wd = 7,ht = 7)


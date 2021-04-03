# --- CrossPath para vias do tipo KO ---

# --- Defina a via, largura e altura 

via <- 'ko00010'
largura <- 1200
altura <- 800
# --- setando o diretorio
setwd("~/CrossPath")
# --- carregando o data proveniente do deseq
#load('data/deseq/gene.rda')

# --- import dos dados da via 
genechoice <- via 
genechoice <- keggGet(genechoice)
genename <- genechoice[[1]][["NAME"]]
# --- transformando o nome da via  
genename <- sapply(strsplit(genename, " / "), "[", 1)
# --- selecionando lista com os genes e seus codigos 
prelista <- as.data.frame(genechoice[[1]][["ORTHOLOGY"]])
prelista <- data.frame(rownames(prelista))
# --- gerando a lista com os codigos dos genes da via
listaNomes <- list()
for (i in 1:nrow(prelista)){
  t <- keggGet(prelista[i,])
  y <- as.data.frame(t[[1]][["GENES"]])
  f <- as.data.frame(y[1,])
  listaNomes <- rbind(listaNomes, f)
}

# --- limpando as linhas que não são HSA(humanos)
listadevias <- as.data.frame(dplyr::filter(listaNomes, grepl("HSA", Genes)))
# --- pegando os genes 
genes_via <- data.frame()
for (i in 1:nrow(listadevias)){
  c <- (ex_between(listadevias[i,], "(", ")"))
  c <- as.data.frame(c[[1]])
  genes_via <- rbind(genes_via, c)
}

# --- comando apenas para deixar 'bonitinho' na lista
colnames(genes_via) <- ("Genes")  

# --- salvar um arquivo com a lista de genes de cada via caso nao tenha sido feito ainda, ignora caso contrario 
myfile<-paste0("all genes of ", genename, ".rda",collapse = NULL)
setwd("~/CrossPath/data/genelist")
if (!file.exists(myfile)) {
  save(genes_via, file = myfile, compress = "xz" )
}

# --- surge como 'deseq'no environment (aproximadamente 1 GB)
# load('data/deseq/gene.rda') 

# --- carrega todas as analises que utilizaremos dos deseq's
deseqnormalizado <- 'GenesNormalizadosDESeq.rda'
setwd("~/CrossPath/data/deseq/")
if (file.exists(deseqnormalizado)) {
  load(deseqnormalizado)
} else {
  # --- carrega todas as analises que utilizaremos dos deseq's
  gene = list()
  gene$vsd =  as.data.frame(assay(deseq$vsd))
  gene$nonNormalized =  as.data.frame(log2(1+counts(deseq$dds, normalized=F)))
  gene$normalized =  as.data.frame(log2(1+counts(deseq$dds, normalized=T)))
  gene$results = as.data.frame(deseq$res)
  
  namedeseq <- paste0("GenesNormalizadosDESeq.rda", collapse = NULL )
  save(gene, file = namedeseq, compress = "xz" )
}

# --- filtrando o dataframe com os genes da via selecionada
agenes <- as.character(genes_via[1:nrow(genes_via),1])
genesFiltered <- list()
genesFiltered$vsd = gene$vsd[unique(agenes),]
genesFiltered$nonNormalized = gene$nonNormalized[unique(agenes),]
genesFiltered$normalized = gene$normalized[unique(agenes),]
genesFiltered$results = gene$results[unique(agenes),]

# --- salvando matrix, caso necessário
genefiltermatrix <- paste0(genename, " matrix", ".rda", collapse = NULL)
setwd("~/CrossPath/data/genematrix")
if (file.exists(genefiltermatrix)){
  load(genefiltermatrix)
} else {
  matrix_result <- genesFiltered$vsd %>% drop_na()
  myfile2<-paste0(genename, " matrix", ".rda", collapse = NULL)
  save(matrix_result, file = myfile2, compress = "xz" )
}

# --- salvando resultados, caso necessário
resultfilter <- paste0(genename, " results", ".rda", collapse = NULL)
setwd("~/CrossPath/data/results")
if (file.exists(resultfilter)){
  load(resultfilter)
} else {
  myfile3<-paste0(genename, " results", ".rda", collapse = NULL)
  save(genesFiltered, file = myfile3, compress = "xz" )
}

# volta para o diretorio inicial
setwd("~/CrossPath")
# --- setando os resultados VSD para montar o grafico
resultadosDaVia <- genesFiltered$results %>% drop_na()  
# --- formatando nomes para criação do gráfico VolcanoPlot
fileName <- paste0("graphics/", genename, " VolcanoPlot", ".png",collapse = NULL)
png(filename = fileName, width = largura, height = altura)
# --- gerando o grafico 
EnhancedVolcano(resultadosDaVia,
                lab = rownames(resultadosDaVia),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                title = paste("VolcanoPlot ", genename, sep = ''))

ggsave(filename = fileName)
graphics.off()
# --- CrossPath para vias do tipo HSA ---

# --- Defina a via, largura e altura 

via <- 'hsa04110'
largura <- 1200
altura <- 800
# --- setando o diretorio
setwd("~/CrossPath")
# --- carregando o data proveniente do deseq
#load('data/deseq/gene.rda')

# --- import dos dados da via 
genechoice <- via  #"hsa04110" 
genechoice <- keggGet(genechoice)
genename <- genechoice[[1]][["NAME"]]
# --- transformando o nome da via  
genename <- sapply(strsplit(genename, " - "), "[", 1)
# --- selecionando lista com os genes e seus codigos 
prelista <- as.data.frame(genechoice[[1]][["GENE"]])
# --- gerando a lista com os codigos dos genes da via
get_cod <- list()
for (i in 1:nrow(prelista)) {
  if (i %% 2){
    nome = data.frame(prelista[i,])
    get_cod = rbind(get_cod, nome)
  }
}
get_cod <- data.frame(get_cod)

# --- adicionando o valor 'hsa' para obter cada gene
cod_final <- list()
j <- "hsa:"
for (i in 1:nrow(get_cod)) {
  cod1<- as.data.frame(get_cod[i,])
  t <- data.frame(lapply(cod1, as.character), stringsAsFactors=FALSE)
  k <- paste(j,t,sep = '')
  k <- as.data.frame(k)
  cod_final <- rbind(cod_final,k) 
  #cod_final contem a lista com todos os codigos para o novo acesso
}

# --- gerando uma lista com todos os possiveis nomes dos genes contidos na via para adequar a qualquer banco de dados
# --- por acessar cada gene, esse loop leva um tempo para capturar todos os dados
lista_all_names <- data.frame()
for (i in 1:nrow(cod_final)) {
  p <- keggGet(cod_final[i,])
  x <- as.data.frame(p[[1]][["NAME"]])
  lista_all_names <- rbind(lista_all_names,x)
}
# --- por fim o resultado é uma lista com todos os nomes de cada gene em uma unica linha
# --- gerando lista com todos os possiveis nomes separados
genes_via <- data.frame()
for (i in 1:nrow(lista_all_names)){
  c <- ex_between(lista_all_names[i,], left = '', right = '(,|$)', fixed = FALSE)
  c <- as.data.frame(c[[1]])
  for (l in 1:nrow(c)){
    if ( l %% 2){
      v <- as.data.frame(c[l,])
      genes_via <- rbind(genes_via, v)
    }
  }
}
# --- comando apenas para deixar 'bonitinho' na lista
colnames(genes_via)<- c("Genes")

# --- salvar um arquivo com a lista de genes de cada via caso nao tenha sido feito ainda, ignora caso contrario 
myfile<-paste0("all genes of ", genename, ".rda",collapse = NULL)
setwd("~/CrossPath/data/genelist")
if (!file.exists(myfile)) {
    save(genes_via, file = myfile, compress = "xz" )
}

# --- surge como 'deseq'no environment (aproximadamente 1 GB)
# load('data/deseq/gene.rda') 

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
####
#
## PCA TFIDF features
##

setwd("~/Documents/KEGG-study/")


df = read.csv("~/Documents/KEGG-study/data/KEGG_querry_tfidf.csv", header = TRUE, stringsAsFactors = FALSE)

df$conf.discrete = "neither"
df$conf.discrete[which(df$over > 0.75)] = "over 0.75"
df$conf.discrete[which(df$over < 0.25)] = "under 0.25"

df.tfidf = df[,grep("pc", names(df))]

pca.df = prcomp(as.matrix(cbind(df.tfidf, 
                                log10(1e-10+apply(df[,grep("HE", names(df))], 1, median)),
                                log10(1e-10+apply(df[,grep("CD.", names(df))], 1, median)),
                                log10(1e-10+apply(df[,grep("UC", names(df))], 1, median)),
                                log10(1e-10+apply(df[,grep("LS", names(df))], 1, median)))))
pca.res = cbind(df$X, df$conf.discrete, as.data.frame(pca.df$x))
names(pca.res)[1] = "kegg"
names(pca.res)[2] = "over.abundance.conf"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = over.abundance.conf)) + geom_point() -> p
ggplotly(p)

ggplot(pca.res, aes(x = PC1, y = PC2, colour = level.1, label = kegg)) + geom_text() -> p
plotly_POST(p, filename = "Keggs-pca-across-keggs", fileopt = "overwrite")
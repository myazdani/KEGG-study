setwd("~/Documents/KEGG-study/")


df = read.csv("./data/keggs_RF_conf_scores.csv", header = TRUE, stringsAsFactors = FALSE)

dist.df = read.csv("./data/tfidf-svd.csv", header = TRUE, stringsAsFactors = FALSE)


first.col.transpose = function(df){
  ## transpose df based on first column as the header
  df.T = as.data.frame(t(as.matrix(df[,-1])))
  names(df.T) = df[,1]
  #df.T$sample.meta = make.names(names(df))
  return(df.T)
}

df.svd = as.data.frame(t(as.matrix(dist.df)))
df.svd$kegg = row.names(df.svd)

df.res = merge(df, df.svd)

library(ggplot2)
library(plotly)


d = dist(as.data.frame(t(as.matrix(dist.df))))
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim


df.mds = cbind(df, as.data.frame(fit$points))

ggplot(df.mds, aes(x = V1, y = V2, colour = log10(1e-8+CD.SRS301867))) + geom_point()


library(tsne)
tsne_data <- tsne(as.data.frame(t(as.matrix(dist.df))), k=2,  max_iter=500, epoch=100)

library(Rtsne)

rtsne_out <- Rtsne(unique(t(as.matrix(dist.df))))

df.tsne = cbind(df[-which(duplicated(t(as.matrix(dist.df)))),], as.data.frame(rtsne_out$Y))

ggplot(df.tsne, aes(x = V1, y = V2, colour = log(over/(1-over)))) + geom_point()



####

numeric.df = df.res[,c(6:68)]
pca.df = prcomp(log10(1e-11 + as.matrix(numeric.df)))
pca.keggs = cbind(df.res$kegg, as.data.frame(pca.df$x))
pca.res = cbind(pca.keggs, df.res$V1)
names(pca.res)[1] = "kegg"
names(pca.res)[length(names(pca.res))] = "V1"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = V1)) + geom_point() -> p





classificaitons = read.csv("./data/classification-results.csv", header = TRUE, stringsAsFactors = FALSE)


bleh = merge(pca.res, classificaitons)
bleh$classifier_output = as.factor(bleh$classifier_output)

ggplot(bleh, aes(x = PC1, y = PC2, colour = log10(1e-8+UC.SRS071970))) + geom_point() -> p



library(Rtsne)
rtsne_out <- Rtsne(unique(numeric.df))
tsne.df = as.data.frame(rtsne_out$Y)
names(tsne.df) = paste0("tsne.", names(tsne.df))
df.tsne = cbind(df.res[-which(duplicated(numeric.df)),], tsne.df)


ggplot(df.tsne, aes(x = tsne.V1, y = tsne.V2, colour = V1)) + geom_point()

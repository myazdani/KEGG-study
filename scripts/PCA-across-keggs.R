####
#
## PCA across keggs
##

setwd("~/Documents/KEGG-IEEE-BigData/")

df = read.csv("./data/table-kegg-clean.csv", header = TRUE, stringsAsFactors = FALSE)


df.kegg.meta = read.csv("./data/kegg_levels_final.csv", header = TRUE, stringsAsFactors = FALSE)

for (i in c(1:nrow(df.kegg.meta))){
  if(df.kegg.meta$level_1[i] == ""){
    df.kegg.meta$level_1[i] = df.kegg.meta$brite_level_1[i]
    df.kegg.meta$level_2[i] = df.kegg.meta$brite_level_2[i]
  }
}

df.meta = merge(df.kegg.meta[,c(1:3,5)], df, all.x = TRUE, all.y = FALSE)
#write.csv(df.meta, file = "./data/kegg-levels-id-subjects.csv", row.names = FALSE, quote = TRUE)
df.meta = read.csv("./data/kegg-levels-id-subjects.csv", header = TRUE, stringsAsFactors = FALSE)


numeric.df = df.meta[,-c(1:5)]
pca.df = prcomp(log10(1e-11 + as.matrix(numeric.df)))
pca.res = cbind(df.meta$X, df.meta$level_1, as.data.frame(pca.df$x))
names(pca.res)[1] = "kegg"
names(pca.res)[2] = "level.1"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = level.1, label = kegg)) + geom_text() -> p
plotly_POST(p, filename = "Keggs-pca-across-keggs", fileopt = "overwrite")

numeric.he = df.meta[,grep("HE", names(df.meta))]
pca.he = prcomp(log10(1e-11 + as.matrix(numeric.he)))
pca.he.res = cbind(df.meta$X, df.meta$level_1, as.data.frame(pca.he$x))

names(pca.he.res)[1] = "kegg"
names(pca.he.res)[2] = "group.1"
ggplot(pca.he.res, aes(x = PC1, y = PC2, colour = group.1, label = kegg)) + geom_text() -> p
ggplot(pca.he.res, aes(x = PC3, y = PC1, colour = group.1, label = kegg)) + geom_point() -> p
plotly_POST(p, filename = "kegg.pca.he", fileopt = "overwrite")

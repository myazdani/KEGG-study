####
#
## PCA across patients
##

setwd("~/Documents/KEGG-study/")

df.patients = read.csv("./data/table-kegg-clean-transpose.csv", header = TRUE, stringsAsFactors = FALSE)
numeric.df = df.patients[,-c(1, ncol(df.patients))]
pca = prcomp(log10(1e-11 + as.matrix(numeric.df)))

pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + geom_text() -> p
ggplotly(p)
plotly_POST(p, filename = "Keggs-pca-across-subjects", fileopt = "overwrite")

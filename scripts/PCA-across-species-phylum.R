####
#
## PCA across species
##

setwd("~/Documents/KEGG-study/")


library(readxl)
df.species = as.data.frame(readxl::read_excel("./data/bin-2013-1015-bae-ref_cov_abundance.xlsx", sheet = 7))

df.species$microbe = paste0("k__", df.species$superkingdom, "p__", df.species$phylum, "c__", df.species$class,
                            "o__", df.species$order, "f__", df.species$family, "g__", df.species$genus, "s__", df.species$species)

names(df.species) = make.names(names(df.species))

# first row is empty
df.species = df.species[-1, ]

# bottom two rows are also empty
df.species = df.species[-c(nrow(df.species)-1, nrow(df.species)),]

sick.patient = "HE.SRS016585"

subjects = names(df.species)[c(grep("LS00", names(df.species)), grep("HE.S", names(df.species)), 
             grep("CD.S", names(df.species)), grep("UC.S", names(df.species)))]

subject.types = subjects
subject.types[grep("LS00", subjects)] = "LS"
subject.types[grep("CD.S", subjects)] = "CD"
subject.types[grep("UC.S", subjects)] = "UC"
subject.types[grep("HE.S", subjects)] = "HE"



numeric.df = as.matrix(df.species[,subjects])


min.val = min(numeric.df[which(numeric.df > 0)])

pca.df = l(log10(.1*min.val + t(numeric.df)))
pca.res = cbind(subjects, subject.types, as.data.frame(pca.df$x))
names(pca.res)[1] = "subjects"
names(pca.res)[2] = "subject.types"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.types, label = subjects)) + geom_text() -> p
#plotly_POST(p, filename = "Species-PCA-keggs-study", fileopt = "overwrite")

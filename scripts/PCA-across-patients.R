####
#
#
## explore reducing KEGGs by their uniqueness (as a way to get rid of a lot of zeros)
##
#
#
##

setwd("~/Documents/KEGG-IEEE-BigData/")

df.patients = read.csv("./data/table-kegg-clean-transpose.csv", header = TRUE, stringsAsFactors = FALSE)
numeric.df = df.patients[,-c(1, ncol(df.patients))]

num.unique = apply(numeric.df, 2, FUN = function(x) length(unique(x)))
numeric.df.reduced = numeric.df[,which(num.unique > 55)]


## naive imputation
#pca = prcomp(log10(1e-9 + as.matrix(numeric.df)))

## prep random imputations:
##

m <- 1e-9
s <- 0
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
draws <- rlnorm(n=63*ncol(numeric.df.reduced), location, shape)
d <- matrix(draws, nrow = 63, byrow = TRUE)
data_imputed = d + as.matrix(numeric.df.reduced)
data_imputed_n = data_imputed/apply(data_imputed, 1, sum)

pca = prcomp(log10(data_imputed_n))



pca.res = cbind(df.patients$subject.type, df.patients$X, as.data.frame(pca$x))
names(pca.res)[2] = "subject.ID"
names(pca.res)[1] = "subject.type"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = subject.type, label = subject.ID)) + 
  geom_point(size = 3) + 
  ggtitle(paste(ncol(data_imputed_n), "KEGGs with imputation mean = ", m, "and SD = ", s)) + 
  theme(legend.title=element_blank())-> p
print(p)


ggplotly(p)

api_create(p, filename = paste0("Keggs_pca_across_subjects_imputation_mean_", m, "_sd_", s), 
           fileopt = "overwrite", sharing = "public")



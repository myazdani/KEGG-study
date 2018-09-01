###

setwd("~/Documents/KEGG-study/")


df = read.csv("./data/keggs_RF_conf_scores.csv", header = TRUE, stringsAsFactors = FALSE)

df.rel = df[,c(1,6:67)]

df.T = as.data.frame(t(as.matrix(df.rel[,-1])))
names(df.T) = df.rel$kegg

df.T$cohort = "HE"
df.T$cohort[grep("LS0", row.names(df.T))] = "LS"
df.T$cohort[grep("CD.", row.names(df.T))] = "CD"
df.T$cohort[grep("UC.", row.names(df.T))] = "UC"
df.T$cohort = as.factor(df.T$cohort)

glm.fit(x = df.T$K00001, y = df.T$cohort, family = "binomial")




aics = list()

for(i in c(1:(ncol(df.T)-1))){
  
  aics[[i]] = glm(cohort ~ ., data = df.T[,c(i,ncol(df.T))], family = "binomial")$aic
  
  #num.sig[[i]] = length(which(x < .01))
}

plot(unlist(aics))
min(unlist(aics))
mean(unlist(aics))

library(plotly)
#[1] 1757 3572 5420 6408  135  441 2822  163 8989 4644  439 7443 5978 8342 8581 8967 4432 3244
#[19] 6945 2896

kegg.df = df.T[,c(135, ncol(df.T))]

ggplot(kegg.df, aes(x = c(1:nrow(kegg.df)), y = log10(1e-6+K00150), colour = cohort)) + geom_point()

ggplot(df.T, aes(x = c(1:nrow(kegg.df)), y = log10(1e-6+K00150), colour = cohort)) + 
  geom_point() +ylab("K00150") -> p

plotly_POST(p, filename = "K00150", fileopt = "overwrite")



ggplot(df.T, aes(x = c(1:nrow(kegg.df)), y = log10(1e-6+K07965), colour = cohort)) + 
  geom_point() +ylab("K07965") -> p

plotly_POST(p, filename = "K07965", fileopt = "overwrite")
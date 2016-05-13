####
#
## kegg analysis
##

setwd("~/Documents/KEGG-study/")

df.meta = read.csv("./data/kegg-levels-id-subjects.csv", header = TRUE, stringsAsFactors = FALSE)


###-------------------------------------------------------------------------------
# KS list
###-------------------------------------------------------------------------------

library(broom)
library(plyr)
df.T = as.data.frame(t(as.matrix(df.meta[,c(6:ncol(df.meta))])))
names(df.T) = df.meta$kegg
df.T$subject = names(df.meta)[c(6:ncol(df.meta))]
df.T$subject.type = "sick"
df.T$subject.type[grep("HE", df.T$subject)] = "healthy"

training.indx = sample(10012,.5*10012)

KEGGs.ks.tests = ldply(lapply(df.T[,-c(training.indx, 10013, 10014)], FUN = function(x) 
  tidy(ks.test(x[which(df.T$subject.type == "healthy")], x[which(df.T$subject.type == "sick")]))))



head(KEGGs.ks.tests[order(KEGGs.ks.tests$statistic, decreasing = TRUE),], 100)


KS.keggs = head(KEGGs.ks.tests[order(KEGGs.ks.tests$statistic, decreasing = TRUE),1], 100)

df.kegg.list = subset(df.meta, kegg %in% KS.keggs)


df.kegg.m = melt(df.kegg.list, id.vars = c(1:5))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"
df.kegg.m$health.state = "healthy"
df.kegg.m$health.state[which(df.kegg.m$subject.type != "HE")] ="sick"

ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("KS KEGGS list")-> p
ggplotly(p)



library(dplyr)

df.kegg.m %>%
  group_by(kegg) %>%
  summarise(median.healthy = median(value[health.state == "healthy"]),
            median.sick = median(value[health.state == "sick"])) %>%
  as.data.frame() -> kegg.abundance

kegg.abundance$healthy.abundance = "under"
kegg.abundance$healthy.abundance[which(kegg.abundance$median.healthy > kegg.abundance$median.sick)] = "over"

print(table(kegg.abundance$healthy.abundance))



training.set = merge(df.kegg.list, kegg.abundance)
library(randomForest)
training.set$healthy.abundance = as.factor(training.set$healthy.abundance)
rf = randomForest(x = training.set[,c(6:68)], y = training.set$healthy.abundance)

rf.preds = predict(rf, df.meta[,-c(1:5)], type = "prob")

df.meta$over = rf.preds[,1]

df.meta$kegg.selected = "1 Used for KS selection"
df.meta$kegg.selected[training.indx] = "0 Test set"
df.meta$kegg.selected[which(df.meta$kegg %in% 
                              kegg.abundance$kegg[which(kegg.abundance$healthy.abundance == "over")])] = "2 over.abundant"

df.meta$kegg.selected[which(df.meta$kegg %in% 
                              kegg.abundance$kegg[which(kegg.abundance$healthy.abundance == "under")])] = "3 under.abundant"

ggplot(df.meta, aes(x = reorder(kegg, -over), y = over, colour =kegg.selected)) + 
  geom_point() + xlab("") + 
  ylab("Confidence of model on healthy having over abundance of kegg")-> p

numeric.df = df.meta[,-c(1:5, 69,70)]
pca.df = prcomp(log10(1e-11 + as.matrix(numeric.df)))
pca.res = cbind(df.meta$X, df.meta$over, as.data.frame(pca.df$x))
names(pca.res)[1] = "kegg"
names(pca.res)[2] = "over.abundant.health.conf"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = over.abundant.health.conf)) + geom_point() -> p
plotly_POST(p, filename = "Keggs-pca-across-keggs", fileopt = "overwrite")

pca.res$over.conv.level = .5
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf > .75)] = "> 0.75"
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf < .25)] = "< 0.25"

ggplot(pca.res, aes(x = PC1, y = PC2, colour = over.conv.level, label = kegg)) + geom_text() -> p
plotly_POST(p, filename = "Keggs-pca-across-keggs-conf", fileopt = "overwrite")



df.kegg.list = subset(df.meta, X %in% c(as.character(head(pca.res$kegg[order(pca.res$over.abundant.health.conf, decreasing = TRUE)], 20)),
                                           as.character(tail(pca.res$kegg[order(pca.res$over.abundant.health.conf, decreasing = TRUE)], 20))))

df.kegg.m = melt(df.kegg.list, id.vars = c(1:5, 69,70))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"

#df.kegg.m$kegg = as.factor(x = df.kegg.m$kegg, levels = KEGGs.RF)
df.kegg.m$abundant = "over"
df.kegg.m$abundant[which(df.kegg.m$X %in% as.character(tail(pca.res$kegg[order(pca.res$over.abundant.health.conf, decreasing = TRUE)], 20)))] = "under"
df.kegg.m$value = as.numeric(df.kegg.m$value)
ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("RF KEGGS from KS moelling") +xlab("") + facet_wrap(~abundant)-> p
ggplotly(p)
plotly_POST(p, filename = "RF KEGGs from KS model", fileopt = "overwrite")

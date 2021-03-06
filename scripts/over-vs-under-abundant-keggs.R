####
#
## kegg analysis
##

setwd("~/Documents/KEGG-study/")

df.meta = read.csv("./data/kegg-levels-id-subjects.csv", header = TRUE, stringsAsFactors = FALSE)
#####
## remove healthy sick subject
#####
df.meta = df.meta[,-ncol(df.meta)]

###-------------------------------------------------------------------------------
# Take random 50% samples and generate KS list
###-------------------------------------------------------------------------------

library(broom)
library(plyr)
df.T = as.data.frame(t(as.matrix(df.meta[,c(6:ncol(df.meta))])))
names(df.T) = df.meta$kegg
df.T$subject = names(df.meta)[c(6:ncol(df.meta))]
df.T$subject.type = "sick"
df.T$subject.type[grep("HE", df.T$subject)] = "healthy"
df.T$cohort = "HE"
df.T$cohort[grep("CD.", df.T$subject)] = "CD"
df.T$cohort[grep("UC.", df.T$subject)] = "UC"
df.T$cohort[grep("LS0", df.T$subject)] = "LS"

set.seed(1000)
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
  ggtitle("KS KEGGS list from 50% random selection")-> p
ggplotly(p)



##---------------------------------------------------------
# label KS Keggs as over or under abundant
##---------------------------------------------------------
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
rf = randomForest(x = training.set[,c(6:67)], y = training.set$healthy.abundance)

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

######### 
# visualize results
#########
training.set.results.from.ks = subset(df.meta, kegg.selected != "0 Test set")
training.over.50 = head(training.set.results.from.ks[order(training.set.results.from.ks$over, decreasing = TRUE),], 50)
training.under.50 = tail(training.set.results.from.ks[order(training.set.results.from.ks$over, decreasing = TRUE),], 50)
training.over.50$group = "TS over"
training.under.50$group = "TS under"
hold.out.set.results = subset(df.meta, kegg.selected %in% c("0 Test set"))
hold.out.set.results = hold.out.set.results[order(hold.out.set.results$over, decreasing = TRUE),]
hold.out.over.50 = head(hold.out.set.results, 50)
hold.out.under.50 = tail(hold.out.set.results, 50)
hold.out.over.50$group = "HS over"
hold.out.under.50$group = "HS under"

rf.res.top50 = rbind(training.over.50, training.under.50, hold.out.over.50, hold.out.under.50)

df.kegg.m = melt(rf.res.top50, id.vars = c(1:5, 68:70))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"

ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("Relative abundance of KEGGs from hold out set and training set") +xlab("") + facet_wrap(~group, ncol = 1)-> p
ggplotly(p)
#plotly_POST(p, filename = "Keggs-from-hold-out-and-train-sets-HE-sick-removed", fileopt = "overwrite")
#plotly_POST(p, filename = "Keggs-from-hold-out-and-train-sets", fileopt = "overwrite")
TS.over = subset(df.kegg.m, group == "TS over")
ggplot(TS.over, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() -> p
ggplotly(p)
TS.under = subset(df.kegg.m, group == "TS under")
ggplot(TS.under,  aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() -> p
ggplotly(p)

HS.over = subset(df.kegg.m, group == "HS over")
ggplot(HS.over,  aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() -> p
ggplotly(p)

HS.under = subset(df.kegg.m, group == "HS under")
ggplot(HS.under,  aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() -> p
ggplotly(p)

bryn.over = c("K03480", "K03483", "K03475", "K02794", "K03753")
bryn.under = c("K01847", "K01711", "K00971", "K12111")


TS.over.list = c("K12111",
                 "K12257",
                 "K13683",
                 "K00971",
                 "K01847",
                 "K00712",
                 "K01206")

HS.over.list = c("K07696",
                 "K01412",
                 "K02199",
                 "K01711",
                 "K02199",
                 "K00957",
                 "K13002")

TS.under.list = c("K03480",
                  "K03483",
                  "K03750",
                  "K03753",
                  "K06351",
                  "K07469",
                  "K09963")

HS.under.list = c("K07757",
                  "K10793",
                  "K02794",
                  "K03475",
                  "K04028",
                  "K06924",
                  "K12527")


# disease.under = rbind(subset(training.set.results.from.ks, kegg %in% TS.over.list), 
#                       subset(hold.out.set.results,  kegg %in% HS.over.list))

melt.helper = function(input.df, col.indx){
  df.kegg.m = melt(input.df, id.vars = col.indx)
  df.kegg.m$variable = as.character(df.kegg.m$variable)
  df.kegg.m$subject.type = "HE"
  df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
  df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
  df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
  df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"
  return(df.kegg.m)
}




disease.under = rbind(subset(df.meta, kegg %in% unique(TS.over.list)[c(1:5)]), 
                      subset(df.meta,  kegg %in% unique(HS.over.list)[c(1:5)]))
disease.under$group = "Hold out set under abundant"
disease.under$group[c(1:5)] = "Training set under abundant" 

disease.under.m = melt.helper(disease.under, c(1:5, 68:70))

ggplot(disease.under.m,  aes(y = kegg, x = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("Disease Under") + facet_wrap(~group, ncol = 1) +ylab("")-> p
ggplotly(p)
#plotly_POST(p, filename = "Sick-patient-under-5-keggs", fileopt = "overwrite")

disease.over = rbind(subset(df.meta, kegg %in% unique(TS.under.list)[c(1:5)]), 
                    subset(df.meta,  kegg %in% unique(HS.under.list)[c(1:5)]))
disease.over$group = "Hold out set over abundant"
disease.over$group[c(1:5)] = "Training set over abundant" 

disease.over.m = melt.helper(disease.over, c(1:5, 68:70))
ggplot(disease.over.m,  aes(y = kegg, x = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("Disease Over") + facet_wrap(~group, ncol = 1) +ylab("")-> p
ggplotly(p)

diseases = rbind(disease.over, disease.under)
diseases.m = melt.helper(diseases, c(1:5, 68:70))
ggplot(diseases.m,  aes(y = kegg, x = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  facet_wrap(~group, ncol = 1) +ylab("")-> p
ggplotly(p)
#plotly_POST(p, filename = "Sick-patient-under-5-keggs", fileopt = "overwrite")

######### 
# PCA style visualization
#########
numeric.df = df.meta[,-c(1:5, 69,70)]
pca.df = prcomp(log10(1e-11 + as.matrix(numeric.df)))
pca.res = cbind(df.meta$X, df.meta$over, as.data.frame(pca.df$x))
names(pca.res)[1] = "kegg"
names(pca.res)[2] = "over.abundant.health.conf"
ggplot(pca.res, aes(x = PC1, y = PC2, colour = over.abundant.health.conf)) + geom_point() -> p
#plotly_POST(p, filename = "Keggs-pca-across-keggs", fileopt = "overwrite")




df.meta$median.HE = 1e-7+ apply(df.meta[,grep("HE.", names(df.meta))], 1, FUN = function(x) mean(x, na.rm = TRUE))
df.meta$median.CD = 1e-7+ apply(df.meta[,grep("CD.", names(df.meta))], 1, FUN = function(x) mean(x, na.rm = TRUE))
df.meta$median.UC = 1e-7+ apply(df.meta[,grep("UC.", names(df.meta))], 1, FUN = function(x) mean(x, na.rm = TRUE))
df.meta$median.LS = 1e-7+ apply(df.meta[,grep("LS0", names(df.meta))], 1, FUN = function(x) mean(x, na.rm = TRUE))

pca.stats = cbind(df.meta$X, df.meta$over, as.data.frame(pca.df$x), df.meta)
names(pca.stats)[1] = "kegg"
names(pca.stats)[2] = "over.abundant.health.conf"
ggplot(pca.stats, aes(x = PC1, y = PC2, colour = log10(median.CD/median.HE))) + geom_point() ->p 
print(p)
ggplot(pca.stats, aes(x = PC1, y = PC2, colour = log10(median.LS/median.HE))) + geom_point() ->p 
print(p)
ggplot(pca.stats, aes(x = PC1, y = PC2, colour = log10(median.UC/median.HE))) + geom_point() ->p 
print(p)

pca.res = pca.stats
pca.res$over.conv.level = .5
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf > .75)] = "> 0.75"
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf < .25)] = "< 0.25"
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf > .98)] = "> 0.98"
pca.res$over.conv.level[which(pca.res$over.abundant.health.conf < .02)] = "< 0.02"

ggplot(pca.res, aes(x = PC1, y = PC2, colour = over.conv.level, label = kegg)) + geom_text() -> p
#plotly_POST(p, filename = "Keggs-pca-across-keggs-conf", fileopt = "overwrite")



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
#plotly_POST(p, filename = "RF KEGGs from KS model", fileopt = "overwrite")

#write.csv(df.meta, file = "./data/keggs_RF_conf_scores.csv", row.names = FALSE, quote = TRUE)





## post pairs


ggplot(df.T, aes(x =  K00648, y = K01993, colour = cohort, label = subject)) + geom_point() -> p
plotly_POST(p, filename = "kegg-pairs-K00648-K01993", fileopt = "overwrite")




min(sapply(df.T[,c(1:10012)], FUN = function(x) min(x[which(x>0)])))

pca.bleh = prcomp(log10(1e-10 + as.matrix(df.T[,c(1:10012)])))
pca.res.bleh = cbind(df.T$cohort, as.data.frame(pca.bleh$x))
names(pca.res.bleh)[1] = "cohort"

ggplot(pca.res.bleh, aes(x = PC1, y = PC2, colour = cohort)) + geom_point() -> p
ggplotly(p)

plotly_POST(p, filename = "Keggs-pca-across-subjects", fileopt = "overwrite")


mediocre.keggs = df.meta[which((df.meta$over > .49) & (df.meta$over < .51)), ]
mediocre.keggs.m = melt(mediocre.keggs[c(1:5), ],  id.vars = c(1:5))
mediocre.keggs.m$variable = as.character(mediocre.keggs.m$variable)
mediocre.keggs.m$value = as.numeric(as.character(mediocre.keggs.m$value))
mediocre.keggs.m$subject.type = "HE"
mediocre.keggs.m$subject.type[grep("LS", mediocre.keggs.m$variable)] = "LS"
mediocre.keggs.m$subject.type[grep("CD.", mediocre.keggs.m$variable)] = "CD"
mediocre.keggs.m$subject.type[grep("UC.", mediocre.keggs.m$variable)] = "UC"
mediocre.keggs.m$subject.type[grep("HE.", mediocre.keggs.m$variable)] = "HE"
ggplot(mediocre.keggs.m, aes(x = log10(1e-8 + value), y = kegg, colour = subject.type)) + geom_point() -> p
plotly_POST(p, filename = "KEGGs with low conf from RF", fileopt = "overwrite")

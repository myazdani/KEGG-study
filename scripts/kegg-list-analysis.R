####
#
## kegg analysis
##

setwd("~/Documents/KEGG-study/")

###-------------------------------------------------------------------------------
# RIGHT FLARE LIST
###-------------------------------------------------------------------------------

kegg.list = read.csv("~/Downloads/keggs_brite_stuff - right_flare_keggs.csv", header = TRUE, stringsAsFactors = FALSE)
df.meta = read.csv("./data/kegg-levels-id-subjects.csv", header = TRUE, stringsAsFactors = FALSE)

KEGGs.right.flare = kegg.list$KEGG.ID
df.kegg.list = subset(df.meta, kegg %in% KEGGs.right.flare)

library(reshape)

df.kegg.m = melt(df.kegg.list, id.vars = c(1:5))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"

ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + 
  geom_point() + ggtitle("Tip of right flare KEGGs") -> p
ggplotly(p)
#plotly_POST(p, filename = "right-flare-keggs", fileopt = "overwrite")

###-------------------------------------------------------------------------------
# 10X list
###-------------------------------------------------------------------------------

library(readxl)

# read_excel reads both xls and xlsx files
df.10Xkeggs = as.data.frame(read_excel("~/Downloads/table-kegg.mehrdad and bryn.xlsx", sheet = 1))
KEGGs.10X = df.10Xkeggs[,1]

df.kegg.list = subset(df.meta, X %in% KEGGs.10X)


df.kegg.m = melt(df.kegg.list, id.vars = c(1:5))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"

ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("10X KEGGS list")-> p
ggplotly(p)
#plotly_POST(p, filename = "10X KEGGs", fileopt = "overwrite")

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

KEGGs.ks.tests = ldply(lapply(df.T[,-c(10013, 10014)], FUN = function(x) 
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

ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("KS KEGGS list")-> p
ggplotly(p)

#plotly_POST(p, filename = "KS KEGGs", fileopt = "overwrite")


###-------------------------------------------------------------------------------
# RF list
###-------------------------------------------------------------------------------
library(randomForest)
df.T$subject.type = as.factor(df.T$subject.type)
rf = randomForest(x = df.T[,-c(10013, 10014)], y = df.T$subject.type, importance = TRUE, ntree = 10000)
importance.df = as.data.frame(rf$importance)
head(importance.df[order(importance.df$MeanDecreaseGini, decreasing = TRUE),])
plot(log10(importance.df$MeanDecreaseGini[order(importance.df$MeanDecreaseGini, decreasing = TRUE)]))
ggplot(df.T, aes(x = log10(1e-10+K02018), y = log10(1e-10+K10900), colour = subject.type)) + geom_point()


KEGGs.RF = head(row.names(importance.df)[order(importance.df$MeanDecreaseGini, decreasing = TRUE)], 100)

df.kegg.list = subset(df.meta, kegg %in% KEGGs.RF)


df.kegg.m = melt(df.kegg.list, id.vars = c(1:5))
df.kegg.m$variable = as.character(df.kegg.m$variable)
df.kegg.m$subject.type = "HE"
df.kegg.m$subject.type[grep("LS", df.kegg.m$variable)] = "LS"
df.kegg.m$subject.type[grep("CD.", df.kegg.m$variable)] = "CD"
df.kegg.m$subject.type[grep("UC.", df.kegg.m$variable)] = "UC"
df.kegg.m$subject.type[grep("HE.", df.kegg.m$variable)] = "HE"

#df.kegg.m$kegg = as.factor(x = df.kegg.m$kegg, levels = KEGGs.RF)
ggplot(df.kegg.m, aes(x = X, y = log10(1e-8+value), colour = subject.type)) + geom_point() + 
  ggtitle("RF KEGGS list") +xlab("")-> p
ggplotly(p)
#plotly_POST(p, filename = "RF KEGGs", fileopt = "overwrite")


######################
#### intersection
######################
intersect(KEGGs.right.flare, KS.keggs)
KEGGs.10X.clean = unname(sapply(KEGGs.10X, FUN = function(x) strsplit(x, split = "\\(")[[1]][1]))
intersect(KEGGs.10X.clean, KS.keggs)

intersect(KEGGs.10X.clean, KEGGs.right.flare)

length(intersect(KEGGs.RF, KS.keggs))

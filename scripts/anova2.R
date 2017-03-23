# Untargeted metabolomics
metadata <- read.csv('../data/mapping.txt', row.names=1,
                     sep='\t')
table <- t(read.table('../data/XCMS.pos.txt'))


table <- table[order(rownames(table)), ]
metadata <- metadata[order(rownames(metadata)), ]
alldata <- merge(table, metadata)
alldata <- merge(table, metadata, by=0, all=TRUE)
rownames(alldata) <- alldata[,"Row.names"]

df = table
potato <- as.factor(alldata$Potato.type)
process <- as.factor(alldata$Processing.type)
pval.df <- data.frame(potato=numeric(0),
                      process=numeric(0),
                      potato.process=numeric(0)))

for(var in colnames(table)){
    k <- as.vector(table[,var])
    res <- anova(lm(k ~ potato + process + potato*process))
    pvals <- res$"Pr(>F)"
}

# Targeted metabolomics
table <- t(read.table('../data/targeted_metabolomics.txt', sep='\t',
                      header=TRUE, row.names=1))

cultivar <- sapply(rownames(table), function(s) unlist(strsplit(s, '-')[1]))
processing <- sapply(rownames(table), function(s) unlist(strsplit(s, '-')[2]))


df = table
potato <- as.factor(alldata$Potato.type)
process <- as.factor(alldata$Processing.type)
pval.df <- data.frame(potato=numeric(0),
                      process=numeric(0),
                      potato.process=numeric(0)))

for(var in colnames(table)){
    k <- as.vector(table[,var])
    res <- anova(lm(k ~ potato + process + potato*process))
    pvals <- res$"Pr(>F)"
}


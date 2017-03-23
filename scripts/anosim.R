library(vegan)

# Untargeted metabolites
metadata <- as.data.frame(read.csv('../data/mapping.txt',
                                   row.names=1,
                                   sep='\t'))
dm <- as.data.frame(read.table('../data/beta_diversity/gower_XCMS.pos.T.txt'))


# @param metamap A data frame of the metadata
# @param field Metadata grouping variable
# @param cat1 Grouping category of interest
# @param cat2 Grouping category of interest
# @param dm Distance Matrix
trim.metadata.pairwise <- function(metamap, field, cat1, cat2, dm){
    idx <- as.logical((metamap[,field]==cat1) + (metamap[,field]==cat2))
    trim.meta <- as.data.frame(metamap[idx,])
    trim.dm <- as.data.frame(dm[idx, idx])
    return(list(trim.meta, trim.dm))
}

# Entire dataset
potato.type.res <- anosim(as.dist(dm), metadata$Potato.type, permutation=9999)
processing.type.res <- anosim(as.dist(dm), metadata$Processing.type, permutation=9999)
color.type.res <- anosim(as.dist(dm), metadata$Potato.Flesh.Color, permutation=9999)


potato.type = c('All Blue', 'All Red', 'Atlantic',
    'Mountain Rose', 'Purple Majesty', 'Russet Burbank')
processing.type = c('Baked', 'Chips', 'French Fries',
    'Microwaved', 'Raw', 'Steamed')

# Potato type anosim
potato.type.df <- data.frame(R=numeric(0),
                             pval=numeric(0))
potato.type.df['All',] <- c(potato.type.res$statistic,
                            potato.type.res$signif)
for(i in seq(1,length(potato.type)-1)){
    for(j in seq(i+1,length(potato.type))){
        res <- trim.metadata.pairwise(metadata,'Potato.type',
                                      potato.type[i],
                                      potato.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        potato.type.res <- anosim(as.dist(trim.dm),
                                  trim.meta$Potato.type,
                                  permutation=9999)
        rowname <- paste(potato.type[i],potato.type[j], sep=', ')
        potato.type.df[rowname,] <- c(potato.type.res$statistic,
                                      potato.type.res$signif)

    }
}
write.table(potato.type.df, '../results/anosim_potato_type.txt', sep='\t')

# Processing type
processing.type.df <- data.frame(R=numeric(0),
                                 pval=numeric(0))
processing.type.df['All',] <- c(processing.type.res$statistic,
                                processing.type.res$signif)
for(i in seq(1,length(processing.type)-1)){
    for(j in seq(i+1,length(processing.type))){
        res <- trim.metadata.pairwise(metadata,'Processing.type',
                                      processing.type[i],
                                      processing.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        processing.type.res <- anosim(as.dist(trim.dm),
                                      trim.meta$Processing.type,
                                      permutation=9999)
        rowname <- paste(processing.type[i],processing.type[j], sep=', ')
        processing.type.df[rowname,] <- c(processing.type.res$statistic,
                                          processing.type.res$signif)

    }
}
write.table(processing.type.df, '../results/anosim_processing_type.txt', sep='\t')


# Flesh color
color.type = c('purple', 'red', 'white')
color.type.df <- data.frame(R=numeric(0),
                            pval=numeric(0))
color.type.df['All',] <- c(color.type.res$statistic,
                           color.type.res$signif)
for(i in seq(1,length(color.type)-1)){
    for(j in seq(i+1,length(color.type))){

        res <- trim.metadata.pairwise(metadata,'Potato.Flesh.Color',
                                      color.type[i],
                                      color.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        color.type.res <- anosim(as.dist(trim.dm),
                                 trim.meta$Potato.Flesh.Color,
                                 permutation=9999)
        rowname <- paste(color.type[i],color.type[j], sep=', ')
        color.type.df[rowname,] <- c(color.type.res$statistic,
                                     color.type.res$signif)
    }
}
write.table(color.type.df, '../results/anosim_flesh_color.txt', sep='\t')

# Targeted metabolomics
dm <- as.data.frame(read.table('../data/beta_diversity/gower_targeted_metabolomics.T.txt'))
metadata <- as.data.frame(read.csv('../data/mapping_full_raw.txt',
                                   row.names=1,
                                   sep='\t'))
potato.type.res <- anosim(as.dist(dm), metadata$Potato, permutation=9999)
processing.type.res <- anosim(as.dist(dm), metadata$Processing, permutation=9999)
color.type.res <- anosim(as.dist(dm), metadata$flesh_color, permutation=9999)


potato.type = c('All Blue', 'All Red', 'Atlantic',
    'Mountain Rose', 'Purple Majesty', 'Russet Burbank')
processing.type = c('Baked', 'Chips', 'French Fries',
    'Microwaved', 'Raw', 'Steamed')

# Potato type anosim
potato.type.df <- data.frame(R=numeric(0),
                             pval=numeric(0))
potato.type.df['All',] <- c(potato.type.res$statistic,
                            potato.type.res$signif)
for(i in seq(1,length(potato.type)-1)){
    for(j in seq(i+1,length(potato.type))){
        res <- trim.metadata.pairwise(metadata,'Potato',
                                      potato.type[i],
                                      potato.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        potato.type.res <- anosim(as.dist(trim.dm),
                                  trim.meta$Potato,
                                  permutation=9999)
        rowname <- paste(potato.type[i],potato.type[j], sep=', ')
        potato.type.df[rowname,] <- c(potato.type.res$statistic,
                                      potato.type.res$signif)

    }
}
write.table(potato.type.df, '../results/anosim_targeted_potato_type.txt', sep='\t')

processing.type = c('Baked', 'Chipped', 'Fried',
                    'Microwaved', 'Raw', 'Steamed')


# Processing type
processing.type.df <- data.frame(R=numeric(0),
                                 pval=numeric(0))
processing.type.df['All',] <- c(processing.type.res$statistic,
                                processing.type.res$signif)
for(i in seq(1,length(processing.type)-1)){
    for(j in seq(i+1,length(processing.type))){

        res <- trim.metadata.pairwise(metadata,'Processing',
                                      processing.type[i],
                                      processing.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        processing.type.res <- anosim(as.dist(trim.dm),
                                      trim.meta$Processing,
                                      permutation=9999)
        rowname <- paste(processing.type[i],processing.type[j], sep=', ')
        processing.type.df[rowname,] <- c(processing.type.res$statistic,
                                          processing.type.res$signif)

    }
}
write.table(processing.type.df, '../results/anosim_targeted_processing_type.txt', sep='\t')


# Potato flesh color
color.type = c('purple', 'red', 'white')
color.type.df <- data.frame(R=numeric(0),
                            pval=numeric(0))
color.type.df['All',] <- c(color.type.res$statistic,
                           color.type.res$signif)
for(i in seq(1,length(color.type)-1)){
    for(j in seq(i+1,length(color.type))){

        res <- trim.metadata.pairwise(metadata,'flesh_color',
                                      color.type[i],
                                      color.type[j],
                                      dm)
        trim.meta <- as.data.frame(res[1])
        trim.dm <- as.data.frame(res[2])
        color.type.res <- anosim(as.dist(trim.dm),
                                      trim.meta$flesh_color,
                                      permutation=9999)
        rowname <- paste(color.type[i],color.type[j], sep=', ')
        color.type.df[rowname,] <- c(color.type.res$statistic,
                                          color.type.res$signif)

    }
}

write.table(color.type.df, '../results/anosim_targeted_flesh_color.txt', sep='\t')

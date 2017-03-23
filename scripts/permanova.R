library(vegan)

# Untargeted analysis
metadata <- as.data.frame(read.csv('../data/mapping.txt',
                                   row.names=1,
                                   sep='\t'))
# dm <- as.data.frame(read.table('../data/beta_diversity/gower_XCMS.pos.T.txt'))

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

adonis(dm ~ Potato.type + Processing.type + Potato.type*Processing.type,
       data=metadata, perm=999)

# Targeted analysis
dm <- as.data.frame(read.table('../data/beta_diversity/gower_targeted_metabolomics.T.txt'))
metadata <- as.data.frame(read.csv('../data/mapping_full_raw.txt',
                                   row.names=1,
                                   sep='\t'))

adonis(dm ~ Potato.type + Processing.type + Potato.type*Processing.type,
       data=metadata, perm=999)

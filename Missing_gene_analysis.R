# nathan dot lazar at gmail dot com

perms <- 1000 # Number of permutations

chrom_dist <- read.table('Chrom_dist_of_losses.txt', header=T,
                         sep='\t', stringsAsFactors=F)
orthos <- read.table('OrhologSet.txt', header=T, sep='\t',
                     stringsAsFactors=F)

# Check for unique gene names in the orthos table
dim(orthos)  # 12193 genes
length(unique(orthos$Associated.Gene.Name)) # 12191 are unique
orthos[duplicated(orthos$Associated.Gene.Name)==T,] # RPP14 and LEKR1 are repeated
orthos[orthos$Associated.Gene.Name=="RPP14",]
# RPP14 has two gene ids but they both start at the same position, the Anole ortholog
# is has two starts that are close to each other. This would come up as a block of two 
# genes, so I removed the second one.
orthos[orthos$Associated.Gene.Name=="LEKR1",]
# LEKR1 has two gene ids with slightly different loci. Again they are in order in the
# gene list and would look like a block, so I removed the second copy.
orthos <- orthos[duplicated(orthos$Associated.Gene.Name)==F,]

# Permutation choosing groups of genes
##################################################################################
n_genes <- sum(chrom_dist$Genes_in_block) # There are 274 genes

chosen_genes <- matrix(NA, nrow=perms, ncol=n_genes)
# allocate empty matrix, genes are stored in columns, permutations in columns

column <- 1 # Start filling at the first column and increment while looping
for (i in 1:nrow(chrom_dist)) {
# Loop over each row in the chrom dist data frame
  line <- chrom_dist[i,]
  h_chr <- line$Human_chr
  block <- line$Genes_in_block
  
  print(paste("Column:", column))

  tofill <- is.na(chosen_genes[,column])
  while(sum(tofill) > 0) {
    # This checks whether the appropriate column of the chosen_genes matrix is 
    # totally full. This will not be full if selection criteria weren't met
    # for one or more of the permutations.

    # Get a 'seed' gene row on the given chromosome for each permutation
    seeds <- sample(which(orthos$Human.Chromosome.Name == h_chr), 
                    sum(tofill), replace=T)
    # Get a block of genes following the seed if the following criteria are met:
    # 1) The genes must be in order on the same human chromosome
    # 2) The orthologous genes in Anole are on the same chromsome (in order)
    # 3) The genes can't have been picked before

    # checking criteria 1 and 2
    check1 <- rep(T, length(seeds)) # all start as true
    check2 <- rep(T, length(seeds))
    if (block != 1) { # No need to check if the block is only 1 gene
      for (j in 1:(block-1)) { # Check that we're on the same chrom for each gene in the block
        check1 <- check1 & 
          orthos$Human.Chromosome.Name[seeds] == orthos$Human.Chromosome.Name[seeds+j]
          # If one of the following genes is on a different chrom, the that element of 
          # check1 will be false
        check2 <- check2 &
          orthos$Anole.Chromosome.Name[seeds] == orthos$Anole.Chromosome.Name[seeds+j]
      }
      check1[is.na(check1)] <- F
      check2[is.na(check2)] <- F      
    }

    # Checking criteria 3. Add in the new gene names and set them to NA if they are
    # repeated
    if (column == 1) { # No need to check criteria 3 if it's the first column
      # just add in gene names when check1 and check2 are true
      for(j in 0:(block-1)) {
        chosen_genes[tofill, column+j][check1 & check2] <- 
          orthos$Associated.Gene.Name[seeds+j][check1 & check2]
      }
    }
    if (column > 1) { 
      # Add new gene names to chosen_genes (if they meet criteria 1 & 2)
      for(j in 0:(block-1)) {
        chosen_genes[tofill, column+j][check1 & check2] <- 
          orthos$Associated.Gene.Name[seeds+j][check1 & check2]
      }

      # Reset to NA if any of the newly added genes are already selected
      reps <- t(apply(chosen_genes, 1, function(x) x[column:(column+block-1)] %in% x[1:(column-1)]))
      if(block > 1) reps <- apply(reps, 1, max)
      chosen_genes[reps==1, column:(column+block-1)] <- NA
    }
  tofill <- is.na(chosen_genes[,column]) # Checking which rows still need to be filled
  print(paste("To Fill:", sum(tofill)))
  }
  column <- column+block
}

# Check for duplicated genes in Rows
apply(chosen_genes, 1, function(x) sum(duplicated(x)))

# OMIM analysis
################################################################################################
# For each permutation (row of chosen_genes) we see how many genes have OMIM ids and descriptions

missing <- read.table('MissingGenes(274).txt', header=F,
                      sep='\t', stringsAsFactors=F)
missing <- missing[,1]

name2id <- read.table('GeneName2OMIMid.txt', header=T, sep='\t',
                      stringsAsFactors=F, na.strings='-')

id2desc <- read.table('OmimID2description.txt', header=T, sep='\t',
                      stringsAsFactors=F, na.strings='', fill=T)
names(id2desc) <- c('OmimID', 'Description')

# Combine two tables matching by the OMIM id
omim <- merge(name2id, id2desc, keep.all=T)

count_omim <- function(genes) {
# counts the number of OMIM ids and descriptions for a set of genes
  ids <- sum(omim$Gene %in% genes & !is.na(omim$OmimID))
  desc <- sum(omim$Gene %in% genes & !is.na(omim$Description))
  return(c(ids, desc))
}

omim_count <- t(apply(chosen_genes, 1, count_omim))
missing_count <- count_omim(missing)
# Two sided p-value calculation
pvalue1 <- (sum(omim_count[,1] > missing_count[1] | 
                  omim_count[,1] < 2*mean(omim_count[,1]) - missing_count[1])) / perms
pvalue2 <- (sum(omim_count[,2] < missing_count[2] | 
               omim_count[,2] > 2*mean(omim_count[,2]) - missing_count[2])) / perms
if(pvalue1 == 0) pvalue1=paste0('<', 1/perms)
if(pvalue2 == 0) pvalue2=paste0('<', 1/perms)

png(file='OMIM_id_dist.png', height=400, width=650)
hist(omim_count[,1], breaks=30, xlim=range(c(missing_count[1], omim_count[,1])),
     main='Distribution of the number of genes with an OMIM disease id',
     xlab='Number of genes out of 274 with OMIM id')
abline(v=missing_count[1], col='red', lwd=2)
legend(197,120, c(paste("Observed count =",missing_count[1]),
                paste("p-value:", pvalue1)), 
       col='red', lty=c(1,0), lwd=2)
abline(v=mean(omim_count[,1]), col='blue', lwd=2)
legend(185,120, paste("Mean count =", mean(omim_count[,2])), 
       col='blue', lty=1, lwd=2, xjust=1)
dev.off()

png(file='OMIM_term_dist.png', height=400, width=650)
hist(omim_count[,2], breaks=30, xlim=range(c(missing_count[2], omim_count[,2])),
     main='Distribution of the number of genes associated \nwith OMIM disease terms',
     xlab='Number of genes out of 274 with OMIM disease term')
abline(v=missing_count[2], col='red', lwd=2)
legend(35,60, c(paste("Observed count =",missing_count[2]),
                 paste("p-value:", pvalue2)), 
       col='red', lty=c(1,0), lwd=2)
abline(v=mean(omim_count[,2]), col='blue', lwd=2)
legend(63,60, paste("Mean count =", mean(omim_count[,2])), 
       col='blue', lty=1, lwd=2)
dev.off()

# Mouse functional screen
###################################################################################
# Get distributions of the number of genes expected to be associated with each
# functional term.

gene2phenoid <- read.table('Gene2PhenoId.txt', header=T, sep='\t',
                           stringsAsFactors=F, fill=T,
                           quote='', comment.char='')

gene2phenoid <- gene2phenoid[,c('GeneName', 'Mouse.Phenotype')]
names(gene2phenoid) <- c('gene', 'phenoId')

phenoid2pheno <- read.table('PhenoId2Pheno.txt', header=T,
                            sep='\t', stringsAsFactors=F, fill=T,
                            quote='', comment.char='')
names(phenoid2pheno) <- c('phenoId', 'pheno_short', 'pheno_long')

# Merge two datasets, only  phenotypes genes associated with them are kept
pheno <- merge(gene2phenoid, phenoid2pheno, all.x=T)
pheno$gene <- toupper(pheno$gene)

# Test the number of genes associated with any phenotype in the missing 
# group vs the permutation groups.

missing_pheno <- pheno[pheno$gene %in% missing,]

# MP:0002169 is "no abnormal phenotype". Besides this one, there are 
# 82 genes associated with a phenotype
missing_pheno_count <- length(unique(missing_pheno[missing_pheno$phenoId!="MP:0002169",]$gene))

#unique(missing_pheno[missing_pheno$phenoId!="MP:0002169",][,1:3])

perm_pheno_count <- rep(0, perms)
for(i in 1:perms) {
  gene_pheno <- pheno[pheno$gene %in% chosen_genes[i,],]
  perm_pheno_count[i] <- length(unique(gene_pheno[gene_pheno$phenoId!="MP:0002169",]$gene))
}
pvalue3 <- sum(perm_pheno_count < missing_pheno_count |
                 perm_pheno_count > 2*mean(perm_pheno_count)-missing_pheno_count)/perms

png(file='Phenotype_count_dist.png', height=400, width=650)
hist(perm_pheno_count, breaks=30, xlim=range(c(missing_pheno_count, perm_pheno_count)),
     main='Distribution of the number of genes with an \nassociated mouse phenotype',
     xlab='Number of genes out of 274 associated with a mouse phenotype')
abline(v=missing_pheno_count, col='red', lwd=2)
legend(95,55, c(paste("Observed count =",missing_pheno_count),
                paste("p-value:", pvalue3)), 
       col='red', lty=c(1,0), lwd=2)
abline(v=mean(perm_pheno_count), col='blue', lwd=2)
legend(133,55, paste("Mean count =", mean(perm_pheno_count)), 
       col='blue', lty=1, lwd=2, xjust=1)
dev.off()

# For each mouse Phenotype with an associated gene (7,874 of them) we find 
# the number of genes in the missing group associated with that phenotype and 
# the number of genes in each permutation associated with that phenotype
# to get a distribution of the expected number of genes associated with each
# phenotype. We then see which phenotypes are significantly represented in the
# missing genes

pheno.df <- data.frame(id = unique(pheno$phenoId),
                       short_desc='', missing=0, mean=0, sd=0, 
                       pvalue=0, long_desc='',
                       stringsAsFactors=F)

for(i in 1:nrow(pheno.df)) {
# Loop through each phenotype id
  id <- pheno.df$id[i]
  sub <- pheno[pheno$phenoId==id,]
  pheno.df$short_desc[i] <- sub$pheno_short[1]
  pheno.df$long_desc[i] <- sub$pheno_long[1]
  
  p_counts <- apply(chosen_genes, 1, function(x) sum(x %in% sub$gene))
  pheno.df$mean[i] <- mean(p_counts)
  pheno.df$sd[i] <- sd(p_counts)
  pheno.df$missing[i] <- sum(missing %in% sub$gene)

  # Calculate p-value depending on which side of the mean the observed 
  # values are on
  if(pheno.df$missing[i] <= pheno.df$mean[i]) {
    pheno.df$pvalue[i] <- sum(p_counts<=pheno.df$missing[i] |
                      p_counts>=(2*pheno.df$mean[i]-pheno.df$missing[i]))/perms
  }
  if(pheno.df$missing[i] > pheno.df$mean[i]) {
    pheno.df$pvalue[i] <- sum(p_counts>=pheno.df$missing[i] |
                                p_counts<=(2*pheno.df$mean[i]-pheno.df$missing[i]))/perms
  }  
}
  
# Sort pheno.df by the number of genes in the missing gene set
# Ties are sorted with the lowest p-value on top
pheno.df <- pheno.df[order(pheno.df$missing, -pheno.df$pvalue, decreasing=T),]

# Save data and Write out table
save(pheno.df.q, file="Mouse_pheno_results.Rdata")  
write.table(pheno.df, file="Mouse_pheno_perm.txt", sep='\t',
            row.names=F)

# Multiple test correction (Benjamini Hochberg)

# No significant genes at the 0.05 q-value level
p <- pheno.df$pvalue
p[p==0] <- 0.001
q <- p.adjust(p, method="fdr")
pheno.df[q<.05,1:6]
sum(q<.05)

# If we assume that values never observed would also not be 
# observed with 10,000 permutations, then we get 47 significant phenotypes
p <- pheno.df$pvalue
p[p==0] <- 0.0001
q <- p.adjust(p, method="fdr")
pheno.df[q<.05,1:6]
sum(q<.05)

# Alternatively we try removing phenotypes where the difference between the mean 
# expected number of genes and the observed number of genes is one or less
pheno.df.q <- pheno.df[abs(pheno.df$mean -pheno.df$missing) > 1,]
nrow(pheno.df.q)
p <- pheno.df.q$pvalue
p[p==0] <- 0.001
q <- p.adjust(p, method="fdr")
pheno.df.q <- cbind(pheno.df.q[,1:6], q, pheno.df.q[,7])
sum(q<.05)

# Save data and Write out table
save(pheno.df.q, file="Mouse_pheno_results_q.Rdata")  
write.table(pheno.df.q, file="Mouse_pheno_perm_q.txt", sep='\t',
            row.names=F)


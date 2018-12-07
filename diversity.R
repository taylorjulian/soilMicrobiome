require(vegan)
require(plyr)
require(lme4)
load('func.Rdata')
load('structure.Rdata')

communityS <- function(x){ # turns abundance tables into community matrices for structure
  list <- mapply(`[<-`, x, 'sample', value = names(x), SIMPLIFY = FALSE)
  bind <- dplyr::bind_rows(list)
  df <- as.data.frame(cast(bind, sample ~ species, value = 'abundance', fun.aggregate = mean))
  df[is.na(df)] <- 0
  df
}

communityF <- function(x){ # turns abundance tables into community matrices for function
  list <- mapply(`[<-`, x, 'sample', value = names(x), SIMPLIFY = FALSE)
  bind <- dplyr::bind_rows(list)
  df <- as.data.frame(cast(bind, sample ~ func, value = 'abundance', fun.aggregate = mean))
  df[is.na(df)] <- 0
  df
}

# alpha diversity
alpha <- function(s, f){
  cs <- communityS(s)
  cf <- communityF(f)
  df <- as.data.frame(cbind(cs[,1], diversity(cs[,-1], index = 'shannon'), diversity(cf[,-1], index = 'shannon')))
  colnames(df) <- c('sample', 'structure', 'func')
  df$structure <- as.numeric(as.character(df$structure))
  df$func <- as.numeric(as.character(df$func))
  df
}

bart <- alpha(structure$bart, func$bart)
cper <- alpha(structure$cper, func$cper)
dsny <- alpha(structure$dsny, func$dsny)
harv <- alpha(structure$harv, func$harv)
jerc <- alpha(structure$jerc, func$jerc)
osbs <- alpha(structure$osbs, func$osbs)
scbi <- alpha(structure$scbi, func$scbi)
ster <- alpha(structure$ster, func$ster)
tall <- alpha(structure$tall, func$tall)
wood <- alpha(structure$wood, func$wood)
alpha <- list(bart, cper, dsny, harv, jerc, osbs, scbi, ster, tall, wood)
names(alpha) <- c('bart', 'cper', 'dsny', 'harv', 'jerc', 'osbs', 'scbi', 'ster', 'tall', 'wood')
alpha <- ldply(alpha)
names(alpha) <- c('site', 'sample', 'structure', 'func')

a <- ggplot(alpha, aes(x = structure, y = func, color = site)) + 
  geom_point() +
  geom_smooth(method = lm)
a1 <- ggplot(alpha, aes(x = structure, y = func)) +
  geom_point() +
  geom_smooth(method = lm)

# beta diversity
beta <- function(s, f, sites){ # function for calculating bray-curtis dissimilarity that places each value in a dataframe
  cs <- communityS(s)
  cf <- communityF(f)
  bs <- as.matrix(vegdist(cs[,-1], 'bray'))
  bf <- as.matrix(vegdist(cf[,-1], 'bray'))
  rownames(bs) <- sites
  colnames(bs) <- sites
  rownames(bf) <- sites
  colnames(bf) <- sites
  dfF <- data.frame(sample1=rownames(bf)[row(bf)], sample2=colnames(bf)[col(bf)],
                        beta.function=c(bf))
  dfS <- data.frame(sample1=rownames(bs)[row(bs)], sample2=colnames(bs)[col(bs)],
                        beta.structure=c(bs))
  df <- merge(dfF, dfS, by = c('sample1', 'sample2'))
  df
}

bart <- beta(structure$bart, func$bart, bartSites)
cper <- beta(structure$cper, func$cper, cperSites)
dsny <- beta(structure$dsny, func$dsny, dsnySites)
harv <- beta(structure$harv, func$harv, harvSites)
jerc <- beta(structure$jerc, func$jerc, jercSites)
osbs <- beta(structure$osbs, func$osbs, osbsSites)
scbi <- beta(structure$scbi, func$scbi, scbiSites)
ster <- beta(structure$ster, func$ster, sterSites)
tall <- beta(structure$tall, func$tall, tallSites)
wood <- beta(structure$wood, func$wood, woodSites)
beta <- list(bart, cper, dsny, harv, jerc, osbs, scbi, ster, tall, wood)
names(beta) <- c('bart', 'cper', 'dsny', 'harv', 'jerc', 'osbs', 'scbi', 'ster', 'tall', 'wood')
beta <- ldply(beta)
names(beta) <- c('site', 'sample1', 'sample2', 'func', 'structure') # dataframe of dissimilarity for all pairwise samples for all sites


b <- ggplot(beta, aes(x = structure, y = func, col = site)) + 
  geom_point() + 
  geom_smooth(method = lm) +
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0, 1) +
  ylim(0, 1)

ggMarginal(b, type = 'density')

# multiple regression with mantel -  are sites close together, functionally similar?
# are site that are functionally similar, structurally similar?
# are sites close together
# three partials, each controlled by the other
# ecodist mhal
# contrast matrix? code each site and code a distance matrix using dist(), and use this as factors for an ANOVA
dist <- distm(samples[,2],samples[,3])


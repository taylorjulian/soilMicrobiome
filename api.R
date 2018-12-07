library(httr)
library(jsonlite)


# function for iterating through different sites structural data
structure <- function(site){
  url <- 'http://api.metagenomics.anl.gov'
  path <- paste0('matrix/organism?group_level=species&source=RefSeq&id=', site)
  response <- GET(url = url, path = path)
  content <- rawToChar(response$content)
  submit <- fromJSON(content)
  response <- GET(submit$url)
  biom <- rawToChar(response$content)
  x <- fromJSON(biom)
  df <- data.frame(x$data$rows$metadata$hierarchy, x$data$data)
  colnames(df)[colnames(df) == 'x.data.data'] <- c('abundance')
  df
}

# function for iterating through different sites functional data
func <- function(site){
  url <- 'http://api.metagenomics.anl.gov'
  path <- paste0('matrix/function?id=', site, '&source=Subsystems&group_level=level1&identity=80')
  response <- GET(url = url, path = path)
  content <- rawToChar(response$content)
  submit <- fromJSON(content)
  response <- GET(submit$url)
  biom <- rawToChar(response$content)
  x <- fromJSON(biom)
  df <- data.frame(x$data$rows$id, x$data$data)
  colnames(df) <- c('func', 'abundance')
  df
}

# vectors of samples
s1 <- c() # of the first type
for(i in 809:868){
  s1[i] <- paste0('mgm4637', i, '.3')
  
}

s2 <- c() # of the second type
for(i in 850:931){
  s2[i] <- paste0('mgm4664', i, '.3')
}

# vectors of samples by site
bartSamples <- c(s1[809])
cperSamples <- c(s1[810:833], s2[858], s2[866], s2[889], s2[890:891], s2[901])
dsnySamples <- c(s1[834:849], s2[865], s2[870], s2[897], s2[910], s2[912], s2[915], 
               s2[918], s2[922], s2[925:926])
jercSamples <- c(s1[850:852])
osbsSamples <- c(s1[853:856], s2[850], s2[852], s2[856:857], s2[861], s2[864], s2[877], 
               s2[882], s2[884], s2[886:887], s2[900], s2[908], s2[913:914], s2[920])
sterSamples <- c(s1[857:868], s2[878], s2[880], s2[898], s2[899], s2[921], s2[924], s2[927], s2[931])
scbiSamples <- c(s2[851], s2[867], s2[896], s2[907])
harvSamples <- c(s2[853:854], s2[862:863], s2[868:869], s2[873], s2[876], s2[881], s2[888], s2[892:893], 
               s2[902:903], s2[909], s2[911], s2[916], s2[919], s2[928])
woodSamples <- c(s2[855], s2[859], s2[871:872], s2[875], s2[883], s2[885], s2[904], s2[917])
tallSamples <- c(s2[860], s2[874], s2[879], s2[894:895], s2[905:906], s2[923], s2[929], s2[930])

## structure ##
# bart
bart <- list()

for(i in bartSamples){
  name <- paste0(i)
  bart[[i]] <- assign(name, structure(i))
}

# cper
cper <- list()

for(i in cperSamples){
  name <- paste0(i)
  cper[[i]] <- assign(name, structure(i))
}

# dsny
dsny <- list()

for(i in dsnySamples){
  name <- paste0(i)
  dsny[[i]] <- assign(name, structure(i))
}

# jerc
jerc <- list()

for(i in jercSamples){
  name <- paste0(i)
  jerc[[i]] <- assign(name, structure(i))
}

# osbs

osbs <- list()

for(i in osbsSamples){
  name <- paste0(i)
  osbs[[i]] <- assign(name, structure(i))
}

# ster
ster <- list()

for(i in sterSamples){
  name <- paste0(i)
  ster[[i]] <- assign(name, structure(i))
}

# scbi
scbi <- list()

for(i in scbiSamples){
  name <- paste0(i)
  scbi[[i]] <- assign(name, structure(i))
}

# harv
harv <- list()

for(i in harvSamples){
  name <- paste0(i)
  harv[[i]] <- assign(name, structure(i))
}

# wood
wood <- list()

for(i in woodSamples){
  name <- paste0(i)
  wood[[i]] <- assign(name, structure(i))
}

 # tall
tall <- list()

for(i in tallSamples){
  name <- paste0(i)
  tall[[i]] <- assign(name, structure(i))
}

# combine structure data into one list
structure <- list(bart, cper, dsny, harv, jerc, osbs, scbi, ster, tall, wood)
names(structure) <- c('bart', 'cper', 'dsny', 'harv', 'jerc', 'osbs', 'scbi', 'ster', 'tall', 'wood')
save(structure, file = 'structure.Rdata')

## function ##
# bart
bart <- list()

for(i in bartSamples){
  name <- paste0(i)
  bart[[i]] <- assign(name, func(i))
}

# cper
cper <- list()

for(i in cperSamples){
  name <- paste0(i)
  cper[[i]] <- assign(name, func(i))
}

# dsny
dsny <- list()

for(i in dsnySamples){
  name <- paste0(i)
  dsny[[i]] <- assign(name, func(i))
}

# jerc
jerc <- list()

for(i in jercSamples){
  name <- paste0(i)
  jerc[[i]] <- assign(name, func(i))
}

# osbs

osbs <- list()

for(i in osbsSamples){
  name <- paste0(i)
  osbs[[i]] <- assign(name, func(i))
}

# ster
ster <- list()

for(i in sterSamples){
  name <- paste0(i)
  ster[[i]] <- assign(name, func(i))
}

# scbi
scbi <- list()

for(i in scbiSamples){
  name <- paste0(i)
  scbi[[i]] <- assign(name, func(i))
}

# harv
harv <- list()

for(i in harvSamples){
  name <- paste0(i)
  harv[[i]] <- assign(name, func(i))
}

# wood
wood <- list()

for(i in woodSamples){
  name <- paste0(i)
  wood[[i]] <- assign(name, func(i))
}

# tall
tall <- list()

for(i in tallSamples){
  name <- paste0(i)
  tall[[i]] <- assign(name, func(i))
}

# combine function data into one list
func <- list(bart, cper, dsny, harv, jerc, osbs, scbi, ster, tall, wood)
names(func) <- c('bart', 'cper', 'dsny', 'harv', 'jerc', 'osbs', 'scbi', 'ster', 'tall', 'wood')
save(func, file = 'func.Rdata')

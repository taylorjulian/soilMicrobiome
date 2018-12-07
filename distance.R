library(PBSmapping)

# distance matrix
  # data frame of samples and lat-lon
  samples <- read.csv('samples.csv')
  
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
  
  # add site name to samples data frame
  sampleSites <- data.frame(c(bartSamples, cperSamples, dsnySamples, jercSamples, osbsSamples, 
                              sterSamples, scbiSamples, harvSamples, woodSamples, tallSamples), 
                            c(rep('bart', length(bartSamples)),
                              rep('cper', length(cperSamples)),
                              rep('dsny', length(dsnySamples)), 
                              rep('jerc', length(jercSamples)), 
                              rep('osbs', length(osbsSamples)),
                              rep('ster', length(sterSamples)), 
                              rep('scbi', length(scbiSamples)),
                              rep('harv', length(harvSamples)),
                              rep('wood', length(woodSamples)),
                              rep('tall', length(tallSamples))))
  colnames(sampleSites) <- c('id', 'site')
  
  samples <- merge(samples, sampleSites, by = 'id')

  # function detecting UTM zone based on longitude
  long2UTMzone <- function(long) {
    (floor((long + 180)/6) %% 60) + 1
  }
  
  # converts lon-lat to UTM
  lonLat2UTM <- function(lon,lat,southern=F){
    
    xy <- cbind(lon,lat)
    
    if(!is.matrix(xy))xy <- matrix(xy,1,2)
    colnames(xy) <- c('X','Y')
    
    zone <- long2UTMzone(xy[,1])
    
    zz <- sort(unique(zone))
    nz <- length(zz)
    
    ll <- xy*0
    
    for(j in zz){
      xyj <- xy[zone == j,]
      if(!is.matrix(xyj))xyj <- matrix(xyj,1)
      colnames(xyj) <- c('X','Y')
      attr(xyj,'projection') <- 'LL'
      
      ll[zone == j,] <- as.matrix(convUL(xyj,km=F,southern=southern))
    }
    ll
  }
  
  
  distance <- function(site){
    df <- samples[samples$site == site,]
    UTM <- lonLat2UTM(df[,3], df[,2])
    df <- data.frame(df$id, UTM)
    colnames(df) <- c('id', 'X', 'Y')
    dist <- as.matrix(dist(df[,2:3]))
    rownames(dist) <- df[,1]
    colnames(dist) <- df[,1]
    dist
  }
  
  mantel.partial(bs, bf, dist)
  
  
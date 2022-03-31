for(file in Sys.glob('/shares/pi_mcmindsr/telomeres/outputs/SRP149363/01_telomerecat/*.csv')) {
  
  telcatout <- read.csv(file)
  read_depth <- as.numeric(read.table(sub('.csv','_nreads.txt',file)))
  species <- XML::xmlToList(XML::xmlParse(file.path('/shares/pi_mcmindsr/raw_data/SRP149363/WGS',sub('.csv','',basename(file)),'metadata.xml')))[[1]][[2]][[2]][[1]]

  if(!exists('alldat')) {
    alldat <- c(telcatout,read_depth=read_depth,species=species)
  } else {
    alldat <- rbind(alldat,c(telcatout,read_depth=read_depth,species=species))
  }

}

write.table(alldat, file='~/SRP149363tels2.txt', sep='\t', quote=F)

i<-0
for(file in Sys.glob('/shares/pi_mcmindsr/telomeres/outputs/SRP224618/01_telomerecat/*.csv')) {
  i <- i+1
  print(i)
  telcatout <- read.csv(file)
  read_depth <- as.numeric(read.table(sub('.csv','_nreads.txt',file)))
  meta <- XML::xmlToList(XML::xmlParse(file.path('/shares/pi_mcmindsr/raw_data/SRP224618/WGA',sub('.csv','',basename(file)),'metadata.xml')))
  species <- meta[[1]][['Attributes']][sapply(meta[[1]][['Attributes']], function(x) 'host' %in% x$.attrs)][[1]][['text']]
  status <- meta[[1]][['Attributes']][sapply(meta[[1]][['Attributes']], function(x) 'Condition' %in% x$.attrs)][[1]][['text']]
    
  if(!exists('alldat')) {
    alldat <- c(telcatout,read_depth=read_depth,species=species,status=status)
  } else {
    alldat <- rbind(alldat,c(telcatout,read_depth=read_depth,species=species,status=status))
  }
  
}

write.table(alldat, file='~/SRP224618tels.txt', sep='\t', quote=F)

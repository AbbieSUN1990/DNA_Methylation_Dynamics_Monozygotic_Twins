###########################################################
################### DMR analysis ##########################
###########################################################
## data:
## The DMR finder in the default setting assumes a R object with two columns: 
## one with Illumina identifiers (Probe_ID) and one with zeros and ones indicating DMP status 
## (Zero: non-DMP and One:DMP)
## To make the DMR finder also available for other types of data one can set 'illumina' to false and use a more generic input:
## Four columns: 
## chr: chromosome (for example "chr1"), 
## pos: genomic position, 
## Status: DMP status (Zero: non-DMP and One:DMP) 
## Probe_ID

## mismatches = maximum number of of allowed non-DMPs within DMRs
## icd = inter CpG distance
## The DMRs can be found in the dmr.allchr object.


DMRfinder = function(data, mismatches = 3, icd = 1000, illumina = TRUE){
  
  #Check for annotation  
  if(illumina == TRUE){
    check  = length(grep("cg", data[, 1])) > 1
    check2  = length(grep("0", data[, 2])) > 1
  
    if(check == FALSE | check2 == FALSE){
      print("Please check your data, is the format right?")
      } else {
        #library(FDb.InfiniumMethylation.hg19)
        library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
        library(dplyr)
        #InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19)
        InfiniumMethylation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
        probesselect = as.data.frame(data, stringsAsFactors = FALSE)
        IM = InfiniumMethylation[match(probesselect[, 1], InfiniumMethylation$Name),]
        probesselect$chr = IM$chr
        probesselect$pos = IM$pos
        }
  } else {
    probesselect = as.data.frame(data)
  }
  
  pb = txtProgressBar(min = 1, max = 22, style = 3)
  
  #DMR finder loop per chromosome
  DMR = function(x, data, mismatches, icd){
    
    setTxtProgressBar(pb,x)
    close(pb)
    
    candidates <- list()
    
    MAXIMUM_REGION_LENGTH = icd # constant, can be adjusted
    
    x = paste("chr", x, sep="")
    chr = probesselect[probesselect$chr==x,]
    order = order(as.numeric(chr$pos))
    chr.sorted = chr[order, ]
    
    last_coordinate = nrow(chr.sorted)
    next_coordinate = 0 # so that a region that has been called will be skipped
    
    for (i in 1:(last_coordinate - 1)) {
      if (i >= next_coordinate) {
        start_location = chr.sorted$pos[i]
        last_visited_status_loc = start_location
        number_of_items = 1
        
        if (chr.sorted$Status[i] == 1) {
          sum_of_ones = 1
          
          # start crawling loop
          for (j in (i+1):last_coordinate) {
            if (chr.sorted$pos[j] > (last_visited_status_loc + MAXIMUM_REGION_LENGTH)) {break}
            if((number_of_items - sum_of_ones) > mismatches) {break}   #Number of mismatches
            number_of_items = number_of_items + 1
            if (chr.sorted$Status[j] == 1) { 
              last_visited_status_loc = chr.sorted$pos[j]
              sum_of_ones = sum_of_ones + 1 
            }
          }
          
          # now check if the result is good enough
          if (sum_of_ones >= 3) {
            last_one = i + number_of_items - 1
            for (k in (i + number_of_items - 1):1) {
              if (chr.sorted$Status[k] == 0) {
                last_one = last_one - 1
                number_of_items = number_of_items - 1
              }
              else {
                break
              }
            }
            
            candidates[[length(candidates) + 1]] = data.frame(
              chr = x,
              start = start_location,
              end = chr.sorted$pos[last_one],
              nDMPs = sum_of_ones,
              N = number_of_items,
              probes = paste(chr.sorted$Probe_ID[chr.sorted$pos >= start_location & chr.sorted$pos <= chr.sorted$pos[last_one] & chr.sorted$Status == 1],
                             collapse = ";"),
              stringsAsFactors = FALSE
            )

            next_coordinate = last_one + 1
            
          }
        }
      }
    }
    
    if(length(candidates) != 0) {
      
      dmr = bind_rows(candidates) %>%
        arrange(start, end) %>%
        distinct()
      
    }
    
    }
  
  #Run DMR finder for all chromosomes
  dmr.func = lapply(seq(1:22), DMR, probesselect, mismatches, icd)
  dmr.allchr = do.call(rbind, dmr.func)
  dmr.allchr
}

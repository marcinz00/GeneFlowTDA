#install.packages("TDA")
#install.packages("fitdistrplus")

library("TDA")
library(fitdistrplus)

setwd("./data")


TDA_MAXDIM <- 1
TDA_MAXSCALE <- 200

#TDA_library <- "PHAT"
TDA_library <- "GUDHI"
#TDA_library <- "Dionysus"

RES_FILE <- "_res_file.csv"



FileProcessing <- function(data_file, res_table)
{
  
  cat("  -> ", data_file, "\n")
  
  # Import raw data:
  data_raw <- readLines(data_file)
  
  
  # find separation line:
  id_sep_0 <- which(data_raw == "BEGIN Network;")
  id_sep_1 <- which(data_raw == "EDGES")
  id_sep_2 <- which(data_raw == "END; [Network]")
  id_sep_3 <- which(data_raw == "VERTICES")
  id_sep_4 <- which(data_raw == "VLABELS")
  
  # create ranges of both data sets:
  #data_1_range <- 4:(id_sep-1)
  #data_2_range <- (id_sep+4):length(data_raw)
  data_0_range <- (id_sep_0 + 1) : (id_sep_0 + 1)
  data_1_range <- (id_sep_1 + 1) : (id_sep_2 - 2)
  data_2_range <- (id_sep_3 + 1) : (id_sep_4 - 2)
  
  # using ranges and row data import it:
  data_0 <- read.table(textConnection(data_raw[data_0_range]), sep = " ")
  data_0[,3] <- gsub("nvertices=", "", data_0[,3])
  NVERT <- as.numeric(data_0[1,3])
  data_0[,4] <- gsub("nedges=", "", data_0[,4])
  data_0[,4] <- gsub(";", "", data_0[,4])
  NEDG <- as.numeric(data_0[1,4])

  
  data_1 <- read.table(textConnection(data_raw[data_1_range]), sep = " ", fill = TRUE, 
                       header = FALSE, col.names=c("c1", "c2", "c3", "c4"), stringsAsFactors = FALSE)
  nrowD1 <- nrow(data_1)
  if ( nrowD1 >= 1 )
    for (i in 1:nrowD1)
    {
      if( nchar(toString((data_1[i,4]))) <= 0 )
        data_1[i,4] <- "w=0.0"
    }
  data_1[,3] <- gsub(",", "", data_1[,3])
  data_1[,4] <- gsub("w=", "", data_1[,4])
  data_1[,4] <- gsub(",", "", data_1[,4])
  
  data_2 <- read.table(textConnection(data_raw[data_2_range]), sep = " ")
  data_2[,3] <- gsub(",", "", data_2[,3])
  

  
  distX <- matrix(.Machine$double.xmax, ncol = NVERT, nrow = NVERT)
  for (i in 1:NVERT)
    distX[i, i] <- 0
  for (i in 1:length(data_1[,1]))
  {
    v1 <- as.numeric( data_1[i, 2] )
    v2 <- as.numeric( data_1[i, 3] )
    
    v1_x <- as.numeric( data_2[v1, 2] )
    v1_y <- as.numeric( data_2[v1, 3] )
    v2_x <- as.numeric( data_2[v2, 2] )
    v2_y <- as.numeric( data_2[v2, 3] )

    dist <- sqrt( (v1_x-v2_x)^2 + (v1_y-v2_y)^2 )
  
    distX[v1, v2] <- dist
    distX[v2, v1] <- dist
  }
  
  
  diag.info <- ripsDiag(X=distX, maxdimension=TDA_MAXDIM, maxscale=TDA_MAXSCALE, 
                        dist="arbitrary", location = TRUE, library=TDA_library, 
                        printProgress=FALSE) 
  
  plot(diag.info$diagram)
  
  ImgFile <- paste(data_file, "_pers.png", sep="")
  dev.copy(png, ImgFile)
  dev.off()
  
  plot(diag.info$diagram, barcode=TRUE)
  
  ImgFile <- paste(data_file, "_bar.png", sep="")
  dev.copy(png, ImgFile)
  dev.off()
  

  h0_ind <- which ( diag.info$diagram[,1] == 0 )
  h0 <- diag.info$diagram[h0_ind,]
  h0_len <- h0[,3] - h0[,2]
  h0_birth <- h0[,2]
  
  h1_ind <- which ( diag.info$diagram[,1] == 1 )
  h1 <- diag.info$diagram[h1_ind,]
  h1_len <- h1[,3] - h1[,2]
  h1_birth <- h1[,2]
  
  h0_mean <- mean(h0_len)
  h1_mean <- mean(h1_len)
  
  h0_median <- median(h0_len)
  h1_median <- median(h1_len)
  
  h0_len_all <- length( h0_len )
  h0_len_25 <- length( which( h0_len > 25 ) )
  
  h1_len_all <- length( h1_len )
  h1_len_100 <- length( which( h1_len > 100 ) )
  h1_birth_100_200 <- length( which( h1_birth >= 100 & h1_birth <=200 ) )
  
  h0_fit.exp <- fitdist(h0_len, "exp")
  plot(h0_fit.exp)

  ImgFile <- paste(data_file, "_bar_exp.png", sep="")
  dev.copy(png, ImgFile)
  dev.off()
  
  h0_exp_lambda <- as.numeric( h0_fit.exp$estimate )
  
    
  row <- c()
  row <- cbind(file_name, h0_mean, h1_mean, h0_mean, h1_median, 
               h0_len_all, h0_len_25, h1_len_all, h1_len_100, h1_birth_100_200, h0_exp_lambda)
  
  res_table <- rbind(res_table, row)

    
  return (res_table)

}


# =========== ALGORITHM =============

result_table <- c()

files <- list.files(path=".", pattern="*.nex", ignore.case = TRUE)
nf <- length(files)
cat("Files:\n")
if ( nf >= 1 )
  for (i in 1:nf)
  {
    file_name <- files[i]
    
    ext <- tolower( substr(file_name, nchar(file_name)-3, nchar(file_name)) )
    if ( ext == ".nex" )
    {
      result_table <- FileProcessing(file_name, result_table)
    }
  }


write.table(result_table, file = RES_FILE, append=FALSE, 
            col.names=TRUE, row.names=FALSE, sep=";")
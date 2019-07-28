## code to prepare simulated data


################################################
#           kink data
################################################
bpoint <- 10000
binwidth <- 50
bins_l <- 40
bins_r <- 40
minz <- bpoint - binwidth*bins_l
maxz <- bpoint + binwidth*bins_r
z <- seq(minz,maxz,binwidth)
bins <- length(z)

maxcount <- 80
mincount <- 30
diff <- (maxcount-mincount)/(bins-1)
count <- seq(maxcount,mincount, -diff)

# put into a dataframe
df <- data.frame(count = count, bin = z)
set.seed(100)
df$count <- df$count + stats::runif(nrow(df),1,15)
df[which(df$bin == bpoint), "count"] <- df[which(df$bin == bpoint), "count"] + 135
res <- lapply(1:nrow(df), function(i) {
    thebin <- df[i,"bin"]
    bincount <- df[i,"count"]
    rep(seq(thebin - (binwidth/2), thebin + (binwidth/2) - 1,10), bincount)
})

kink_vector <- unlist(res)

length(kink_vector)

################################################
#           notch data
################################################

bpoint = 10000
binwidth = 50
bins_l = 40
bins_r = 40
minz <- 8000
maxz <- 12000
z <- seq(minz,maxz,binwidth)
bins <- length(z)

maxcount <- 70
mincount <- 40
diff <- (maxcount-mincount)/(bins-1)
count <- seq(maxcount,mincount, -diff)


# put into a dataframe
df <- data.frame(count = count, bin = z)
set.seed(100)
df$count <- df$count + stats::runif(81,1,23)

df[which(df$bin == bpoint), "count"] <- df[which(df$bin == bpoint), "count"] + 335

# create notch
i <- 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 40
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 30
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 20
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 20
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 20
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 20

i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 20

i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 4
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 4
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 3

i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 3

i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 2
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 2
i <- i + 1
df[which(df$bin == bpoint + i*binwidth), "count"] <- df[which(df$bin == bpoint + i*binwidth), "count"] - 2



df[which(df$bin == 11100), "count"] <- df[which(df$bin == 11100), "count"] - 10

df[which(df$bin == 11150), "count"] <- df[which(df$bin == 11100), "count"] - 10
df[which(df$bin == 11050), "count"] <- df[which(df$bin == 11100), "count"] - 10
df[which(df$bin == 10900), "count"] <- df[which(df$bin == 11100), "count"] - 10
df[which(df$bin == 10700), "count"] <- df[which(df$bin == 11100), "count"] - 10


res <- lapply(1:nrow(df), function(i) {
    thebin <- df[i,"bin"]
    bincount <- df[i,"count"]
    rep(seq(thebin - (binwidth/2), thebin + (binwidth/2) - 1,10), bincount)
})

notch_vector <- unlist(res)

# put into a dataframe of equal length

obs_diff <- length(notch_vector) - length(kink_vector)


bunching_data <- data.frame(kink_vector = c(kink_vector,rep(NA,obs_diff)),
                            notch_vector = notch_vector)


usethis::use_data(bunching_data, overwrite = TRUE, compress = "bzip2")

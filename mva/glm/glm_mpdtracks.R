# Read a txt file, named "variables.txt"

mpdtracks <- read.table("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/variables.txt", 
                      sep ="", header = FALSE, dec =".")

mpdtracks02 <- read.table("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/variables02.txt", 
                        sep ="", header = FALSE, dec =".")

mpdtracks03 <- read.table("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/variables03.txt", 
                          sep ="", header = FALSE, dec =".")


mpdtracks04 <- read.table("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/variables04.txt", 
                          sep ="", header = FALSE, dec =".")

mpdtracks05 <- read.table("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/variables05.txt", 
                          sep ="", header = FALSE, dec =".")


result<-file("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/mpdtracks_results.txt")

result<-file("/home/julio/Research/2018-2021/00_must/MexNICA/analysis/R/mpdtracks_results02.txt")

names(mpdtracks)[1] <- "eta"
names(mpdtracks)[2] <- "P"
names(mpdtracks)[3] <- "dEdx"
names(mpdtracks)[4] <- "PDGID"

names(mpdtracks02)[1] <- "eta"
names(mpdtracks02)[2] <- "P"
names(mpdtracks02)[3] <- "dEdx"
names(mpdtracks02)[4] <- "PDGID"

names(mpdtracks03)[1] <- "eta"
names(mpdtracks03)[2] <- "P"
names(mpdtracks03)[3] <- "dEdx"
names(mpdtracks03)[4] <- "PDGID"

names(mpdtracks04)[1] <- "eta"
names(mpdtracks04)[2] <- "P"
names(mpdtracks04)[3] <- "dEdx"
names(mpdtracks04)[4] <- "PDGID"

names(mpdtracks05)[1] <- "eta"
names(mpdtracks05)[2] <- "P"
names(mpdtracks05)[3] <- "dEdx"
names(mpdtracks05)[4] <- "PDGID"

head(mpdtracks)

class(mpdtracks)

str(mpdtracks)


model <- glm(formula = PDGID ~ dEdx + P, data=mpdtracks, family = binomial)

summary(model)

model

pm=NULL
ym=NULL
z=NULL

j = 0




for (i in 1:20000){
  newdata = data.frame(P = mpdtracks05$P[i], dEdx = mpdtracks05$dEdx[i]) 
  x = predict(model, newdata, type="response")
  pm[i] = x
  if(x <= 0.5){
    ym[i] = 0
    if( (ym[i]-mpdtracks05$PDGID[i]) == 0){
      z[i] = 100
    }
    else z[i] = 0
  }
  else if(x > 0.5){
    ym[i] = 1
    if( (ym[i]-mpdtracks05$PDGID[i]) == 0){
      z[i] = 100
    }
    else z[i] = 0
  }
  if(z[i]==0) j = j + 1
}

#df <- data.frame(
 # mpdtracks$PDGID,
  #ym, 
  #z
#)

print(j)

#write.table(df, result, append = FALSE, sep = " ", dec = ".",
 #           row.names = FALSE, col.names = TRUE)

close(result)






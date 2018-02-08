## SCRIPT TO PREP DATA FOR INTERMEDIATE R WORKSHOP

#install tidyverse

#install.packages("tidyverse",repos = "http://ftp.ussg.iu.edu/CRAN/")
install.packages("nycflights13", repos = "http://ftp.ussg.iu.edu/CRAN/")

library(tidyverse)
library(nycflights13)

capwords <- function(s, strict = TRUE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                  {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


housing<-read.csv("https://raw.githubusercontent.com/datasets/house-prices-us/master/data/cities-month.csv", stringsAsFactors=F, strip.white = T)
housing=housing[c(1:(length(housing)-3),length(housing))]
airpass<-as.data.frame(t(matrix(AirPassengers,12,dimnames=list(month.abb,unique(floor(time(AirPassengers)))))))
airpass$Year = rownames(airpass)
rownames(airpass)<-NULL
airpass<-airpass[c(13,1:12)]
nzmodeshareschool<-read.csv("http://www.transport.govt.nz/assets/CSVFiles/TP007Modeshareofjourneystoschoolbyregionage51220102014.csv", stringsAsFactors=F)
nzmodeshareschool <- nzmodeshareschool[1:8,1:10]
names(nzmodeshareschool)[1] = c("travel_mode")
mms<-read.table("http://www.randomservices.org/random/data/MM.txt", header=TRUE,stringsAsFactors=F)
mms$BagID = seq(1:length(mms$Red))
mms <- mms %>% tbl_df %>% gather(color, count, -Weight, -BagID) %>% arrange(BagID,color)
canlang<-read.csv("http://www12.statcan.ca/census-recensement/2011/dp-pd/hlt-fst/lang/Tables/Download/FileDownload.cfm?TabID=1&Lang=E&Asc=1&PRCode=01&OrderBy=1&View=1&tableID=401&queryID=2&Age=1&DLType=2&Delim=1", skip=2, strip.white = T, stringsAsFactors = F)
aircodes<-read.csv("https://raw.githubusercontent.com/datasets/airport-codes/master/data/airport-codes.csv", stringsAsFactors = F)
aircodes <- tbl_df(aircodes)
us_aircodes <- aircodes %>% filter(iso_country %in% c("US", "PR", "VI", "AR"), type != "closed", iata_code != "")
uncity<-read.csv("https://raw.githubusercontent.com/datasets/population-city/master/data/unsd-citypopulation-year-fm.csv", stringsAsFactors=F)
uscity <- uncity %>% tbl_df %>% filter(Country.or.Area == "United States of America") %>% select(City, Year, Sex, Value) %>% group_by(City, Year) %>% summarize(pop = sum(as.numeric(Value))) %>% separate(City, c("city", "state"), sep="\\(") %>% mutate(state = trimws(sub(")", "", state)), city=apply(as.matrix(trimws(city),ncol=1),1,capwords), Year=as.numeric(Year))

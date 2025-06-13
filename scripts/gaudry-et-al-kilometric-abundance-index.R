
setwd("C:/Users/jeany/OneDrive/Documents/ofb")
ab1 <-read.table("data_abondance.txt", sep="\t", header =T, dec=".", na.strings = "")

df.all <- data.frame()

for(i in 2:ncol(ab1)){

coef.tp <- coef(lm(ab1[,i]~ab1$annee))[2]

confint.tp <- confint(lm(ab1[,i]~ab1$annee))[2,]

all <- c(coef.tp, confint.tp)

df.all <- rbind(df.all, all)
}

df.all <- df.all * 100
df.all <- round(df.all,2)
colnames(df.all) <- c("trend","ic.low","ic.up")
write.csv2(df.all,"percent_change.csv")

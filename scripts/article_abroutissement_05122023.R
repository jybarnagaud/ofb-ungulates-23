#--------------------------------------------------#
####### R script for Gaudry et al. [TITLE] #########
# contact : william.gaudry@ofb.gouv.fr
# last modified : 5/12/2023
# replicates the analyses presented in the paper and 
# performs some additional technical checks.
#--------------------------------------------------#

### !!! raw data files to be limited to necessary variables
# and column labels to be translated to english before releasing this script

### for figures : homogenize names of variables and order of presentation of variables
# for coherence with the text. 

### for figures : !!! homogenize scales within figures!!!


#---------------------------#
#### necessary libraries ####
#---------------------------#

library(ade4)
library(PerformanceAnalytics)
library(tidyr)
library(ggplot2)
library(lme4)
library(mgcv)
library(ggeffects)
library(patchwork)
library(Hmisc)
library(adegraphics)
library(gratia)
library(viridis)
library(cowplot)
library(DHARMa)
library(mapview)
library(sf)
library(spdep)
library(mgcViz)
library(sjPlot)
library(oddsratio)

#------------#
#### Data ####
#------------#

# all data (extracted by Mathieu Garel - OFB : september 2023). Browsing is 
# 1 (browsed) or 0 (not browsed) for each plot - year

load("data/pour_jy_will.RData")
abr <- data

# data for species-level analyses (extracted by Colombe Lefort - OFB : july 2023)
# browsing is for plot - year - species

app.sp <- read.csv2("data/scores_appetence.csv")
div.pla <- read.csv2("data/Diversite_vegetale_placette.csv")

## prepare spatial data  -------------------------------------------------------

spat.data <- unique(
  abr[, c(
    "Numero.placette",
    "Massif",
    "Densite_lin",
    "Frequentation",
    "Altitude",
    "Pente",
    "Tirs",
    "Distance",
    "dist.trail",
    "dist.road",
    "dist.linear",
    "length.linear",
    "length.road",
    "length.trail",
    "strava.buff",
    "slope",
    "elev",
    "northness",
    "vrm",
    "slope.buff",
    "vrm.buff",
    "northness.buff",
    "elev.buff",
    "vis.buff"
  )])

## prepare temporal and spatial-temporal data ----------------------------------

spatemp.data <-
  abr[, c("ID_unique",
          "Massif",
          "Annee",
          "Numero.placette",
          "Appetence_mean",
          "Div_spe_veg","Bois")]

## Browsing data : remove missing plant species names --------------------------

div.pla1 <- subset(div.pla,Nom_Latin!="Aucune")

abr.data0 <-
  unique(div.pla1[, c("Nom_Latin",
                     "Annee",
                     "Numero.placette",
                     "Presence",
                     "Consommation")])

## Browsing data : aggregate data ----------------------------------------------
# in some instances a given species is shown as browsed and non-browsed on a plot
# a given year. If so, it is considered as browsed (same for presence). This has
# been cross checked with field operators. 

abr.data1 <- aggregate(abr.data0[,c("Presence","Consommation")],by = list(abr.data0$Nom_Latin,
                                                    abr.data0$Annee,
                                                    abr.data0$Numero.placette),FUN= "sum")

abr.data1[which(abr.data1$Presence >= 1),"Presence"] = 1

colnames(abr.data1) = c("Nom_Latin","Annee","Numero.placette","Presence","Consommation")

abr.data <- subset(abr.data1,Numero.placette %in% unique(abr$Numero.placette))
abr.data$index <- paste(abr.data$Numero.placette,abr.data$Annee,sep="_")

#-----------------------------#
#### Prepare presence data ####
#-----------------------------#

## Presence matrix (to retrieve absences) --------------------------------------

pres.long <- abr.data[,c("Presence","index","Nom_Latin")]
pres.wide <- as.data.frame(pivot_wider(pres.long,id_cols = index,values_from = Presence,names_from = Nom_Latin))
pres.wide[is.na(pres.wide)] <- 0
rownames(pres.wide) <- pres.wide$index
pres.wide1 <- pres.wide[,-1]
plots.id <- factor(substring(rownames(pres.wide1),first=1,last=5))

rows <- rowSums(pres.wide1)
cols <- colSums(pres.wide1)

# remove species absent from all plots

abr.data2 <- subset(abr.data,Nom_Latin %nin% names(which(cols==0)) & Annee %in% c(2007:2014))

pres.long <- abr.data2[,c("Presence","index","Nom_Latin")]
pres.wide <- as.data.frame(pivot_wider(pres.long,id_cols = index,values_from = Presence,names_from = Nom_Latin))
pres.wide[is.na(pres.wide)] <- 0
rownames(pres.wide) <- pres.wide$index
pres.wide1 <- pres.wide[,-1]
plots.id <- factor(substring(rownames(pres.wide1),first=1,last=5))

rows <- rowSums(pres.wide1)
cols <- colSums(pres.wide1)

## Browsing matrix -------------------------------------------------------------

abr.long <- abr.data2[,c("Consommation","index","Nom_Latin")]
abr.wide <- as.data.frame(pivot_wider(abr.long,id_cols = index,values_from = Consommation,names_from = Nom_Latin))
abr.wide[is.na(abr.wide)] <- 0
rownames(abr.wide) <- abr.wide$index
abr.wide1 <- abr.wide[,-1]

rows <- rowSums(abr.wide1)
cols <- colSums(abr.wide1)

# Browsing is considered as binary

abr.sum <- apply(abr.wide1,1,sum)
abr2 <- abr.sum
abr2[abr2>0] <- 1

#--------------------------#
#### Explore covariates ####
#--------------------------#

## correlations among spatial variables ----------------------------------------

spat.data2 <- as.data.frame(
  spat.data[, c(
    "Numero.placette",
    "Massif",
    "elev.buff",
    "slope.buff",
    "northness.buff",
    "vrm.buff",
    "dist.linear",
    "strava.buff",
    "Tirs",
    "vis.buff"
  )]
)

spat.data3 <- spat.data2
rownames(spat.data3) <- spat.data3$Numero.placette
spat.data3 <- spat.data3[unique(as.character(plots.id)),]

chart.Correlation(spat.data3[,-c(1:2,ncol(spat.data3))])

## PCA on the subset of variables kept -----------------------------------------

pca.spat.data3 <- dudi.pca(spat.data3[,-c(1,2,ncol(spat.data3))],scannf=F,nf=3)

screeplot(pca.spat.data3)
s.corcircle(pca.spat.data3$co)
s.corcircle(pca.spat.data3$co,xax = 1, yax = 3)
s.label(pca.spat.data3$li)
s.class(pca.spat.data3$li,fac=factor(spat.data3$Massif))

#-----------------------------#
#### Community level model ####
#-----------------------------#

## Prepare community-level data ------------------------------------------------

# compute the species richness of present and browsed species

abr.rs <- apply(abr.wide1,1,sum)
pres.rs <- apply(pres.wide1,1,sum)

# retrieve plot number and year

Numero.placette <- substring(names(abr.rs), first = 1, last = 5)
Annee <-
  as.numeric(substring(names(abr.rs), first = 7, last = nchar(names(abr.rs))))
df.vr <-
  data.frame(
    Numero.placette = Numero.placette,
    Annee = Annee,
    conso = abr.rs,
    presence = pres.rs
  )

# merge covariate, presence and browsing data

df.vrs <-
  merge(df.vr, spat.data3[, c(
    "Numero.placette",
    "Massif",
    "elev.buff",
    "northness.buff",
    "vrm.buff",
    "dist.linear",
    "Tirs",
    "vis.buff"
  )], by = "Numero.placette", all = F)

# join with temporal variables - reduces the datasets to the years actually sampled
# for each plot (n = 2584 plots x year, 411 plots and 8 years)

df.vrst0 <-
  merge(df.vrs,
        spatemp.data[, c("Numero.placette", "Annee", "Appetence_mean", "Bois")],
        by = c("Numero.placette", "Annee"),
        all = F)

df.vrst <- as.data.frame(df.vrst0)
df.vrst$Massif <- factor(df.vrst$Massif)
df.vrst$Bois <- factor(df.vrst$Bois)

## explanatory variables  ------------------------------------------------------

# check distributions

par(mfrow=c(4,3))
hist(df.vrst$elev.buff)
hist(df.vrst$northness.buff)
hist(df.vrst$vrm.buff)
hist(df.vrst$dist.linear)
hist(df.vrst$Tirs)
hist(df.vrst$Appetence_mean)
barplot(table(df.vrst$Bois))
hist(log(df.vrst$dist.linear+1))
hist(log(df.vrst$Tirs+1))
hist((df.vrst$vis.buff))

# regularization (log transform)

df.vrst$Ldist.linear <- log(df.vrst$dist.linear+1)
df.vrst$LTirs <- log(df.vrst$Tirs+1)
df.vrst$Lpresence <- log(df.vrst$presence)

# stand type as a factor

df.vrst$Bois <- factor(df.vrst$Bois)

# browsing as a binary variable

df.vrst$conso.bin <- df.vrst$conso
df.vrst[which(df.vrst$conso.bin > 0),"conso.bin"] <- 1

## Descriptive figure of browsing pressure -------------------------------------

descr.dat <-
  aggregate(df.vrst$conso.bin,
            by = list(df.vrst$Annee, df.vrst$Massif),
            FUN = "sum")
eff.dat <-
  aggregate(
    df.vrst$Numero.placette,
    by = list(df.vrst$Annee, df.vrst$Massif),
    FUN = "length"
  )

all.descr <- merge(eff.dat, descr.dat, by = c("Group.1", "Group.2"))

colnames(all.descr) <- c("year", "massif", "n.plot", "n.brows")

all.descr$prop.brows = 100 * all.descr$n.brows / all.descr$n.plot

# do the plot : variation per year

supp.labs <- c("(b) Hautes Bauges", "(c) Semnoz", "(d) Cimeteret")
names(supp.labs) <- c("HAUTES BAUGES", "SEMNOZ", "SW BAUGES")

descr.curves <- ggplot(all.descr) +
  aes(x = year, y = prop.brows) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ massif, labeller = labeller(massif = supp.labs)) +
  theme_classic() +
  theme(
    strip.background = element_blank() ,
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    )
  ) +
  labs(x = "Years", y = "Proportion of browsed plots (%)") 
# boxplot : variation between regions

descr.regions <- ggplot(all.descr) +
  aes(x = massif, y = prop.brows) +
  geom_boxplot(fill = "gray90") +
  theme_classic() +
  labs(x = "Cluster", y = "Proportion of browsed plots (%)") +
  scale_x_discrete(
    labels = c(
      "HAUTES BAUGES" = "Hautes Bauges",
      'SEMNOZ' = 'Semnoz',
      'SW BAUGES' = 'Cimeteret'
    )
  ) +
  theme(
    strip.background = element_blank() ,
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("(a) Regional variation in browsing pressure")
    
# both

combined_plot <- descr.regions / descr.curves 
combined_plot

ggsave("outputs/curves_prop_brows.png", width = 7, height = 10 )

## Generalized Additive Model --------------------------------------------------

k <- 4

abrbin.glob.gam <-
  gam(
    conso.bin ~ s(elev.buff, k = k) + 
      s(northness.buff, k = k) +
      s(vrm.buff, k = k) + 
      s(Ldist.linear, k = k) + 
      s(LTirs, k = k) + 
      s(Appetence_mean, k = k) + 
      s(vis.buff,k = k)+
      s(Annee, k = k) +
      s(Lpresence, k = k)+
     Bois + 
      Massif ,
    family = binomial,
    data = df.vrst
  )

# check GAM residuals 
gam.check(abrbin.glob.gam)
summary(abrbin.glob.gam)

## Generalized Linear Model on scaled variables --------------------------------

abrbin.glob.lm <-
  glm(
    conso.bin ~ scale(elev.buff) + 
      scale(northness.buff) +
      scale(vrm.buff) + 
      scale(Ldist.linear) + 
      scale(LTirs) + 
      scale(Appetence_mean) + 
      scale(vis.buff) + 
      scale(Annee) +
      scale(Lpresence)+
      Bois + 
      Massif ,
    family = binomial,
    data = df.vrst
  )

## GLM on non-scaled variables -------------------------------------------------

abrbin.nosc.glob.lm <-
  glm(
    conso.bin ~ elev.buff + 
      northness.buff +
      vrm.buff + 
      Ldist.linear + 
      LTirs + 
      Appetence_mean + 
      vis.buff + 
      Annee +
      Lpresence +
      Bois + 
      Massif ,
    family = binomial,
    data = df.vrst
  )

# store coefficients

sc.coefs.lm <- names(coef(abrbin.glob.lm)[2:10])

# dispersion diagnostics  

sim.gam <-
  simulateResiduals(fittedModel = abrbin.glob.gam, plot = F)
testDispersion(sim.gam)

sim.glm.sc <-
  simulateResiduals(fittedModel = abrbin.glob.lm, plot = F)
testDispersion(sim.gam)

sim.glm <-
  simulateResiduals(fittedModel = abrbin.nosc.glob.lm, plot = F)
testDispersion(sim.gam)

# spatial autocorrelation diagnostic (on the GAM)

xy <- unique(div.pla[,c("Numero.placette","Long","Lat")])

resid.gam <-
  aggregate(residuals(abrbin.glob.gam),
            by = list(df.vrst$Numero.placette),
            FUN = "mean")

resid.gam.xy <-
  merge(resid.gam,
        xy,
        by.x = "Group.1",
        by.y = "Numero.placette",
        all = F)

colnames(resid.gam.xy)[c(1, 2)] <-
  c("Numero.placette", "residuals.gam")

mapView(
  resid.gam.xy,
  xcol = "Long",
  ycol = "Lat",
  zcol = "residuals.gam",
  crs = "epsg:27572"
)

gam.xy <- as.matrix(resid.gam.xy[, c("Long", "Lat")])
dn.gam <- dnearneigh(gam.xy, d1 = 0, d2 = 5000)
corr.gam <-
  sp.correlogram(
    dn.gam,
    var = resid.gam.xy$residuals.gam,
    method = "I",
    order = 4,
    zero.policy = T
  )
plot(corr.gam)

# temporal autocorrelation diagnostic

resid.gam.yr <-
  aggregate(residuals(abrbin.glob.gam),
            by = list(df.vrst$Annee),
            FUN = "mean")
colnames(resid.gam.yr) <- c("Annee", "resid.mean.yr")
acf(resid.gam.yr$resid.mean.yr, main = "residual temporal autocorrelation")

## Table of odds ratios (non-scaled model) -------------------------------------

# biologically relevant increments

increments <-
  list(
    elev.buff = 100,
    northness.buff =  0.25,
    vrm.buff = 0.0005,
    Ldist.linear = log(100),
    LTirs = log(4),
    Appetence_mean = 0.5,
    vis.buff = mean(df.vrst$vis.buff),
    Annee = 1,
    Lpresence = log(3)
  )

# odds ratios

ci.all0 <-
  or_glm(data = df.vrst, model = abrbin.nosc.glob.lm, incr = increments)

ci.all <- as.data.frame(ci.all0[, c(2, 3, 4)])
colnames(ci.all) <- c("odds", "lower", "upper")

ci.all$species <- "community"
ci.all$variable <-
  c(
    "elev",
    "northness",
    "vrm",
    "ldist.lin",
    "ltirs",
    "app_mean",
    "lvisi",
    "annee",
    "sr",
    "GB",
    "PB",
    "Semnoz",
    "SW Bauges"
  )

# Odds ratio table for export

write.table(ci.all,"outputs/table_odds_ratios_community.txt",sep=";")

## forest plots of coefficients with the scaled model [FIGURE 2]----------------

sc.coefs.lm <- names(coef(abrbin.glob.lm)[2:10])

plot_model(
  abrbin.glob.lm,
  terms = sc.coefs.lm,
  colors = "viridis",
  sort.est = T,
  transform = NULL
) +
  theme_classic() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "goldenrod3") +
  labs(x = "Scaled variables", y = "Scaled estimates - 95%CI", title = "")

 ggsave("outputs/fit_scaled_global.png",height=8,width=4)

 
 ## Partial residuals plots of GAM models = [FIGURE 3]--------------------------
 
 p1 <- draw(abrbin.glob.gam, residuals = T, select = 1) +
   labs(x = "Elevation (m)", "Marginal effect", title = "") +
   theme_classic()
 
 p2 <- draw(abrbin.glob.gam, residuals = T, select = 2) +
   labs(x = "Northness (relative index)", "Marginal effect", title = "") +
   theme_classic()
 
 p3 <- draw(abrbin.glob.gam, residuals = T, select = 3) +
   labs(x = "Rugosity (relative index)", "Marginal effect", title = "") +
   theme_classic()
 
 p4 <- draw(abrbin.glob.gam, residuals = T, select = 4) +
   labs(x = "Distance to the nearest linear element (m, log-transformed)", "Marginal effect", title =
          "") +
   theme_classic()
 
 p5 <- draw(abrbin.glob.gam, residuals = T, select = 5) +
   labs(x = "Number of hunting shots (log-transformed)", "Marginal effect", title =
          "")+
     theme_classic()
 
 p6 <- draw(abrbin.glob.gam, residuals = T, select = 6) +
   labs(x = "mean appetency", "Marginal effect", title = "") +
   theme_classic()
 
 p7 <- draw(abrbin.glob.gam, residuals = T, select = 7) +
   labs(x = "Visibility (unit???)", "Marginal effect", title = "") +
   theme_classic()
 
 p8 <- draw(abrbin.glob.gam, residuals = T, select = 8) +
   labs(x = "Years", "Marginal effect", title = "") +
   theme_classic()
 
 p9 <- draw(abrbin.glob.gam, residuals = T, select = 9) +
   labs(x = "Species richness (log-transformed)", "Marginal effect", title =
          "") +
   theme_classic()
 
 par.terms <- parametric_effects(abrbin.glob.gam)
 
 cowplot::plot_grid(
   p1,
   p2,
   p3,
   p4,
   p5,
   p6,
   p7,
   p8,
   p9,
   nrow = 3,
   ncol = 3,
   labels = "AUTO"
 )
 
 ggsave("outputs/fit_gam_global.png",
        height = 15,
        width = 15)

#------------------------------------------#
#### Multi-species comparative analysis ####
#------------------------------------------#
 
## prepare data ----------------------------------------------------------------

# select the 10 most common species (highest number of occurrences)
 
rank.sp <- rev(sort(tapply(abr.data2$Presence,INDEX = abr.data2$Nom_Latin,FUN = "sum")))

sp.keep <- names(rank.sp[1:10])

abr.data3 <- subset(abr.data2,Nom_Latin %in% sp.keep)

# proportion of browsed sites per species

prop.abr <- aggregate(abr.data3[,c("Presence","Consommation")],by=list(abr.data3$Nom_Latin,abr.data3$Annee),FUN="sum")
colnames(prop.abr)[c(1,2)] <- c("Nom_Latin","Annee")

prop.abr$freq.abr <- 100 * (prop.abr$Consommation / prop.abr$Presence)

# graphical displays 

p1 <- ggplot(prop.abr)+
  aes(x = Nom_Latin,  y = freq.abr)+
  geom_boxplot()+
  theme_classic()+
  labs(x = "Species", y = "% grazed plots per year")+ 
  theme(axis.text.x = element_text(angle = 90))

p2 <- ggplot(prop.abr)+
  aes(x = Annee, y = freq.abr)+
  geom_line()+
  theme_classic()+
  labs(x = "Years", y = "% grazed plots")+
  facet_wrap(~Nom_Latin)

p1+p2

ggsave("outputs/explo_multisp.png",width = 10, height = 5)

## multi-species mixed GAM (not used in the article) ---------------------------

abr.data4 <-
  merge(abr.data3,
        df.vrst,
        by = c("Numero.placette", "Annee"),
        all = F)

abr.data5 <-
  merge(abr.data4,
        app.sp,
        by.x = "Nom_Latin",
        by.y = "espece",
        all = F)

abr.data5$Nom_Latin <- factor(abr.data5$Nom_Latin)
abr.data5$Numero.placette <- factor(abr.data5$Numero.placette)

#abrbin.sp.gamm <-
  "gamm(
    Consommation ~ s(elev.buff, k = k, by = Nom_Latin) + 
      s(northness.buff, k = k, by = Nom_Latin) +
      s(vrm.buff, k = k, by = Nom_Latin) + 
      s(Ldist.linear, k = k, by = Nom_Latin) + 
      s(LTirs, k = k, by = Nom_Latin) + 
      s(Appetence_mean, k = k, by = Nom_Latin) + 
      s(Annee, k = k, by = Nom_Latin) +
      s(presence, k = k, by = Nom_Latin)+
      Bois +
      (Bois : Nom_Latin) +
      Nom_Latin +
      Massif +
     (Massif : Nom_Latin),
    random = list(Numero.placette = ~1 ),
    family = binomial,
    data = abr.data5
  )"

#save(abrbin.sp.gamm,"abrbin.sp.gamm.RData")

## species per species GAM -----------------------------------------------------

# Faster to run, avoids a random effects and also easier to represent graphically

# this loop runs models and store results

for(i in 1:nlevels(abr.data5$Nom_Latin)){

sp <- levels(abr.data5$Nom_Latin)[i]

sub.sp <- subset(abr.data5,Nom_Latin == sp)

nm.sp <- paste("conso_gam",sub(" ", "_", sp),sep="_")

# Run one GAM per species 

assign(nm.sp,
 gam(
  Consommation ~ s(elev.buff, k = k) + 
    s(northness.buff, k = k) +
    s(vrm.buff, k = k) + 
    s(Ldist.linear, k = k) + 
    s(LTirs, k = k) + 
    s(Appetence_mean, k = k) + 
    s(vis.buff, k = k)+
    s(Annee, k = k) +
    s(Lpresence, k = k)+
    Bois +
    Massif, 
  family = binomial,
  data = sub.sp
)
)

# run the corresponding GLM with scaled variables 

nm2.sp <- paste("conso_glm",sub(" ", "_", sp),sep="_")
nm2.nosc.sp <- paste("conso_nosc_glm",sub(" ", "_", sp),sep="_")

assign(nm2.sp,
       glm(
         Consommation ~ scale(elev.buff) + 
           scale(northness.buff) +
           scale(vrm.buff) + 
           scale(Ldist.linear) + 
           scale(LTirs) + 
           scale(Appetence_mean) +
           scale(vis.buff) +
           scale(Annee) +
           scale(Lpresence) +
           Bois +
           Massif, 
         family = binomial,
         data = sub.sp
       )
)

# run the corresponding GLM with non-scaled variables (for odds ratios)

assign(nm2.nosc.sp,
       glm(
         Consommation ~ elev.buff + 
           northness.buff +
           vrm.buff + 
           Ldist.linear + 
           LTirs + 
           Appetence_mean + 
           vis.buff + 
           Annee +
           Lpresence +
           Bois +
           Massif, 
         family = binomial,
         data = sub.sp
       )
)
}

# store results in lists

gam.sp.list <- ls()[grep("conso_gam",ls())]
glm.sp.list <- ls()[grep("conso_glm",ls())]
glm.nosc.sp.list <- ls()[grep("conso_nosc_glm",ls())]

## comparison of species responses to each variable [FIGURE 4] -----------------

# compare the curves(uses gratia)

comp.sp.all <-
  compare_smooths(
    get(gam.sp.list[1]),
    get(gam.sp.list[2]),
    get(gam.sp.list[3]),
    get(gam.sp.list[4]),
    get(gam.sp.list[5]),
    get(gam.sp.list[6]),
    get(gam.sp.list[7]),
    get(gam.sp.list[8]),
    get(gam.sp.list[9]),
    get(gam.sp.list[10]),
    smooths = c(
      "s(elev.buff)",
      "s(northness.buff)",
      "s(vrm.buff)",
      "s(Ldist.linear)",
      "s(LTirs)",
      "s(Appetence_mean)",
      "s(vis.buff)",
      "s(Annee)",
      "s(Lpresence)"
    )
  )

# graphical displays (note : there are warnings about NAs which are only due to
# the way the data objects are restructured to perform the plots)

comp.sp.all.df <- unnest(comp.sp.all,cols="data")

p.gam.ele <- ggplot(comp.sp.all.df) +
  aes(x = elev.buff, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Scaled marginal effect", title =
         "") 

p.gam.north <- ggplot(comp.sp.all.df) +
  aes(x = northness.buff, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,3) + 
  labs(x = "Northness (relative unit)", y = "Scaled marginal effect", title =
         "")
  
p.gam.vrm <- ggplot(comp.sp.all.df) +
  aes(x = vrm.buff, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,2) +
labs(x = "Rugosity (relative unit)", y = "Scaled marginal effect", title =
         "")
  
p.gam.lin <- ggplot(comp.sp.all.df) +
  aes(x = Ldist.linear, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
    labs(x = "Distance to nearest linear element (log(m))", y = "Scaled marginal effect", title =
         "")

p.gam.tirs <- ggplot(comp.sp.all.df) +
  aes(x = LTirs, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "number of hunting shots (log)", y = "Scaled marginal effect", title =
         "")

p.gam.app <- ggplot(comp.sp.all.df) +
  aes(x = Appetence_mean, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "mean Appetence", y = "Scaled marginal effect", title =
         "") 

p.gam.vis <- ggplot(comp.sp.all.df) +
  aes(x = vis.buff, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "Visbility (Unit???)", y = "Scaled marginal effect", title =
         "") 

p.gam.an <- ggplot(comp.sp.all.df) +
  aes(x = Annee, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-2, 2) +
  labs(x = "Years", y = "Scaled marginal effect", title =
         "")

p.gam.sr <- ggplot(comp.sp.all.df) +
  aes(x = Lpresence, y = est, col = model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1, 3) +
  labs(x = "Species richness", y = "Scaled marginal effect", title =
         "")

## figure ----------------------------------------------------------------------

cowplot::plot_grid(
  p.gam.ele,
  p.gam.north,
  p.gam.vrm,
  p.gam.lin,
  p.gam.tirs,
  p.gam.app,
  p.gam.vis,
  p.gam.an,
  p.gam.sr,
  nrow = 3,
  ncol = 3,
  labels = "AUTO"
)

ggsave("outputs/species_level_plots.png", width = 9, height = 9)

## compute odds ratios ---------------------------------------------------------

# loops over species to compute odds ratios on non-scaled GLMs

or.glm <- data.frame()

for (i in 1:length(glm.nosc.sp.list)) {
  sp.mod <- get(glm.nosc.sp.list[i])
  
# odds ratios

dat.or <-
  subset(abr.data5, Nom_Latin == levels(abr.data5$Nom_Latin)[i])
elev.sp.or <-
  or_glm(data = dat.or, model = sp.mod, incr = increments)

# cr?er un jeu de donn?es

ci.sp <- as.data.frame(elev.sp.or[, c(2, 3, 4)])
colnames(ci.sp) <- c("odds", "lower", "upper")
ci.sp$species <-
  sub("_", " ", substring(glm.sp.list[i], 11, nchar(glm.sp.list[i])))
ci.sp$variable <-
  c(
    "elev",
    "northness",
    "vrm",
    "ldist.lin",
    "ltirs",
    "app_mean",
    "visi",
    "annee",
    "sr",
    "GB",
    "PB",
    "Semnoz",
    "SW Bauges"
  )
or.glm <- rbind(or.glm, ci.sp)
}

# concatenate with the odds ratios of the community model

ci.sp.all <- rbind(or.glm,ci.all)

# plots odds ratios as a forest plot

varinames <- data.frame(names = c("elev","northness","vrm","ldist.lin","ltirs","app_mean","annee","sr"),
varilabs = c("Elevation (m)","Northness (relative unit","Rugosity (relative unit)","Distance to nearest linear element (log(m))","Number of hunting shots (log)","Mean attractivity (scores)","Years","Species richness")
)

for(i in 1:nrow(varinames)){

vari <- varinames[i,1]
ci.sp.elev <- subset(ci.sp.all,variable == vari & species!="community")
ci.sp.elev$species <- reorder(ci.sp.elev$species,ci.sp.elev$odds)
ci.sp.elev$n_species <- as.numeric(ci.sp.elev$species)

ci.comm.elev <- subset(ci.sp.all,variable == vari & species=="community")
ci.comm.elev$n_species <- 0

assign(paste("plot",varinames[i,1],sep="_"), ggplot(data=ci.sp.elev, aes(y=n_species, x=odds, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  geom_point(aes(x=odds,y=0),data = ci.comm.elev,colour="goldenrod3")+
  geom_errorbarh(aes(xmin=lower,xmax=upper),data = ci.comm.elev,colour="goldenrod3",height=.1)+
  theme_classic()+
  scale_y_continuous(breaks = 0:10,labels = c("community",as.character(levels(ci.sp.elev$species))))+
  labs(x = "Odds ratio", y = "Species",title = varinames[i,2])+
  geom_vline(xintercept = 1, linetype="dashed",col="steelblue")
)
}

# plot of "stand" effect (contrasts)

or.bois.sp <- subset(ci.sp.all, variable%in%c("GB","PB")& species!="community")
or.bois.sp$species <- reorder(or.bois.sp$species,or.bois.sp$odds)
or.bois.sp$n_species <- as.numeric(or.bois.sp$species)

or.bois.comm <- subset(ci.sp.all, variable%in%c("GB","PB")& species=="community")
or.bois.comm$n_species <- 0

plot_bois <- ggplot(data=or.bois.sp, aes(y=n_species, x=odds, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  geom_point(aes(x=odds,y=0),data = or.bois.comm,colour="goldenrod3")+
  geom_errorbarh(aes(xmin=lower,xmax=upper),data = or.bois.comm,colour="goldenrod3",height=.1)+
  theme_classic()+
  scale_y_continuous(breaks = 0:10,labels = c("community",as.character(levels(ci.sp.elev$species))))+
  labs(x = "Odds ratio (by contrast with MB)", y = "Species",title = "Stand type")+
  geom_vline(xintercept = 1, linetype="dashed",col="steelblue") +
    facet_wrap(~variable)
  
# plot of "region" effect (contrasts)

or.zone.sp <- subset(ci.sp.all, variable%in%c("Semnoz","SW Bauges")& species!="community")
or.zone.sp$species <- reorder(or.zone.sp$species,or.zone.sp$odds)
or.zone.sp$n_species <- as.numeric(or.zone.sp$species)

or.zone.comm <- subset(ci.sp.all, variable%in%c("Semnoz","SW Bauges")& species=="community")
or.zone.comm$n_species <- 0

plot_zone <- ggplot(data=or.zone.sp, aes(y=n_species, x=odds, xmin=lower, xmax=upper)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  geom_point(aes(x=odds,y=0),data = or.zone.comm,colour="goldenrod3")+
  geom_errorbarh(aes(xmin=lower,xmax=upper),data = or.zone.comm,colour="goldenrod3",height=.1)+
  theme_classic()+
  scale_y_continuous(breaks = 0:10,labels = c("community",as.character(levels(ci.sp.elev$species))))+
  labs(x = "Odds ratio (by contrast with Hautes Bauges)", y = "Species",title = "Stand type")+
  geom_vline(xintercept = 1, linetype="dashed",col="steelblue") +
  facet_wrap(~variable)

## forest plot of scaled effects [FIGURE 2]-------------------------------------

# has to be combined with the first part of figure 2 above

names(glm.sp.list) <- levels(abr.data5$Nom_Latin)

for (i in 1:length(glm.sp.list)){
  
glm.tp <- get(glm.sp.list[i])
sc.coefs.lm <- names(coef(glm.tp)[2:10])

assign(paste("sp_forplot",i,sep="_"),plot_model(
  glm.tp,
  terms = sc.coefs.lm,
  colors = "viridis",
  sort.est = T,
  transform = NULL
) +
  theme_classic() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "goldenrod3") +
  labs(x = "Scaled variables", y = "Scaled estimates - 95%CI", title = "")
)
}


cowplot::plot_grid(
  sp_forplot_1,
  sp_forplot_2,
  sp_forplot_3,
  sp_forplot_4,
  sp_forplot_5,
  sp_forplot_6,
  sp_forplot_7,
  sp_forplot_8,
  sp_forplot_9,
  sp_forplot_10,
  nrow = 3,
  ncol = 3,
  labels = "AUTO"
)

ggsave("outputs/species_level_forest_plots.png", width = 9, height = 9)



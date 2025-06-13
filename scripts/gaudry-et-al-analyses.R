#------------------------------------------------------------------------------#
            ####### R script for Gaudry et al. [TITLE]  #########
# main author : William Gaudry - Office Français de la Biodiversité
# script author : Jean-Yves Barnagaud - Ecole Pratique des Hautes Etudes
# contact for article and data : william.gaudry@ofb.gouv.fr
# contact for script : jean-yves.barnagaud@ephe.psl.eu
# last modified : 13/06/2025
# replicates the analyses presented in the paper and 
# performs some additional technical checks.
#------------------------------------------------------------------------------#

## libraries -------------------------------------------------------------------


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
#library(devtools)
#install_github("gavinsimpson/gratia",force = T) # needed for some graphical formatting details
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
library(colorspace) 
library(modelsummary)
library(openxlsx)
library(tmap)

## DATA ========================================================================

# all data (extracted by Mathieu Garel - OFB : september 2023). Browsing is 
# 1 (browsed) or 0 (not browsed) for each plot - year

load("data/pour_jy_will.Rdata") 
abr <- data

# data for species-level analyses (extracted by Colombe Lefort - OFB : july 2023)
# browsing is for plot - year - species

app.sp <- read.csv2("data/scores_appetence.csv")
div.pla <- read.csv2("data/Diversite_vegetale_placette.csv")

# Names of variables (full)

variable.names <- c("Species richness (log-transformed)", 
  "Mean appetency",
  "Years",
  "Elevation (m)",
  "Number of hunting shots (log-transformed)",
  "Distance to the nearest linear element (m, log-transformed)",
  "Rugosity (relative index)",
  "Visibility (unit???)",
  "Northness (relative index)")

## PREPARE DATA ================================================================

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

## EXPLORE COVARIATES ==========================================================

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

# note : strava is removed due to its very strong triangular relation with dist.linear

pca.spat.data3 <- dudi.pca(spat.data3[,-c(1,2,8,ncol(spat.data3))],scannf=F,nf=3)

screeplot(pca.spat.data3)
s.corcircle(pca.spat.data3$co)
s.corcircle(pca.spat.data3$co,xax = 1, yax = 3)
s.label(pca.spat.data3$li)
s.class(pca.spat.data3$li,fac=factor(spat.data3$Massif))

## COMMUNITY-LEVEL MODEL =======================================================

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
df.vrst$Massif <- factor(df.vrst$Massif,levels = c("SEMNOZ","SW BAUGES","HAUTES BAUGES"))
df.vrst$Bois <- factor(df.vrst$Bois,levels = c("PB","BM","GB"))

## Figure 1 : export for qGis mapping ------------------------------------------

# the maps from fig. 1 itself are compiled in qGIS. geographic data sources:

# extent of forests in fig. c-d-e : https://geoservices.ign.fr/bdforet (accessed 08/08/2024)
# Bauges Regional Nature Park : https://www.data.gouv.fr/fr/datasets/parcs-naturels-regionaux-de-france-au-15-septembre-2021-58-pnr/ (accessed 08/08/2024)
# cities : https://www.data.gouv.fr/fr/datasets/villes-de-france/ (accessed 08/08/2024)
# the shapefile for the extent of continental France is a standard low resolution, unvalidated polygon

# compute number of sampling years / plot and time range

dur.plot <-
  aggregate(df.vrst$Annee,
            by = list(df.vrst$Numero.placette),
            FUN = "length")

yrange.plot <-
  aggregate(df.vrst$Annee,
            by = list(df.vrst$Numero.placette),
            FUN = range)

yrange.plot<-data.frame(yrange.plot$Group.1,unlist(yrange.plot$x))


ytime <- merge(dur.plot,yrange.plot, by = 1)
colnames(ytime) <- c("Numero.placette","duration","min.year","max.year")

ymap <- merge(ytime,spat.data3,by = "Numero.placette",all = F)
ymap.sf <- st_as_sf(ymap)
# st_write(ymap.sf, dsn = "outputs/fig1_map.geojson", layer = "nc.geojson")

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
hist((df.vrst$presence))

# regularization (log transform)

df.vrst$Ldist.linear <- log(df.vrst$dist.linear+1)
df.vrst$LTirs <- log(df.vrst$Tirs+1)
df.vrst$Lpresence <- log(df.vrst$presence)
df.vrst$Lvrm.buff <- log(df.vrst$vrm.buff)

# stand type as a factor

df.vrst$Bois <- factor(df.vrst$Bois)

# browsing as a binary variable

df.vrst$conso.bin <- df.vrst$conso
df.vrst[which(df.vrst$conso.bin > 0),"conso.bin"] <- 1

## Figure 2a : descriptive figure of browsing pressure --------------------------

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


all.descr$massif <- factor(all.descr$massif,levels = c("SEMNOZ","SW BAUGES","HAUTES BAUGES"))

# boxplot : variation between regions

f2a <- ggplot(all.descr) +
  aes(x = massif, y = prop.brows) +
  geom_boxplot(aes(fill = massif),alpha = 0.4) +
  scale_fill_manual(values = c('#46327e', '#1fa187', 'goldenrod'),name = "cluster", labels = c("Semnoz","Cimeteret","Hautes Bauges"))+
  ylim(0,100)+
  theme_classic() +
  labs(x = "", y = "browsed plots (% per year)") +
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
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill="none")

# curves : variations within regions

f2b <- ggplot(all.descr) +
  aes(x = year, y = prop.brows, group = massif) +
  geom_line(aes(color = massif)) +
  geom_point(aes(color = massif)) +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(0, 100) +
  theme_classic() +
  theme(
    strip.background = element_blank() ,
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 12,hjust = 0)
  ) +
  labs(title = "a",x = "", y = "browsed plots (% per year)")

# all
#combined_plot <- (f2a / f2b)+
 # plot_annotation(tag_levels = "a",
  #                theme = theme(plot.title = element_text(size = 12, face = "bold")))

#ggsave("outputs/box_curves_prop_brows.svg",width = 1750,height = 2500,units = "px")

# only curves (current version in article)

#f2b
#ggsave("outputs/curves_prop_brows.svg",width = 1500,height = 1000,units = "px")


## Figure 2b : descriptive trends per species -----------------------------------

# add community-level browsing to the data table

com.conso <-
  as.data.frame(df.vrst[, c("Numero.placette", "Annee", "conso.bin")])

abr.data2b <- merge(abr.data2,com.conso,by = c("Numero.placette","Annee"))

# select the 10 most common species (highest number of occurrences)

rank.sp <- rev(sort(tapply(abr.data2b$Presence,INDEX = abr.data2b$Nom_Latin,FUN = "sum")))

sp.keep <- names(rank.sp[1:10])

# keep only the 10 most common species and plots with at least 1 grazed species

abr.data3 <- subset(abr.data2b,Nom_Latin %in% sp.keep & conso.bin == 1)

# proportion of browsed sites per species

prop.abr <- aggregate(abr.data3[,c("Presence","Consommation")],by=list(abr.data3$Nom_Latin,abr.data3$Annee),FUN="sum")
colnames(prop.abr)[c(1,2)] <- c("Nom_Latin","Annee")

prop.abr$freq.abr <- 100 * (prop.abr$Consommation / prop.abr$Presence)

# median proportion of browsed plots per species

med.per.sp <- tapply(prop.abr$freq.abr,INDEX = prop.abr$Nom_Latin, FUN = "median")
sort(med.per.sp)

# distributions of yearly browsing pressures per species (not displayed in article)

p1 <- ggplot(prop.abr)+
  aes(x = reorder(Nom_Latin,freq.abr,median),  y = freq.abr)+
  geom_boxplot(fill = "gray90", alpha = 0.5)+
  theme_classic()+
  labs(x = "Species", y = "% browsed plots per year")+ 
  theme(axis.text.x = element_text(face = "italic", angle = 90))

# curves of yearly browsing pressures per species 

prop.abr$Nom_Latin <- reorder(prop.abr$Nom_Latin,prop.abr$freq.abr,median)
labs.sp.yr <- paste(letters[2:(nlevels(prop.abr$Nom_Latin)+1)],levels(prop.abr$Nom_Latin),sep = " - ")
names(labs.sp.yr) <- levels(prop.abr$Nom_Latin)


p2b <- ggplot(prop.abr)+
  aes(x = Annee, y = freq.abr)+
  geom_line()+
  theme_classic()+
  labs(x = "", y = "% browsed plots")+
  facet_wrap(~Nom_Latin,nrow = 4,scales = "free",labeller = as_labeller(labs.sp.yr))+
  theme_classic() +
  theme(
    strip.background = element_blank() ,
    strip.text = element_text(face = "bold", size = 12,hjust=0)
    )+
  ylim(0,100)

f2b + p2b

layout <- "
AAA
AAA
AAA
CCC
CCC
CCC
CCC
"
(f2b / p2b) + 
  plot_layout(design = layout)

ggsave(
  "outputs/gaudry_et_al_fig2.svg",
  width = 2500,
  height = 4000,
  units = "px"
)

## Generalized Additive Model --------------------------------------------------

k <- 4

abrbin.glob.gam <-
  gam(
    conso.bin ~ s(elev.buff, k = k) + 
      s(northness.buff, k = k) +
      s(Lvrm.buff, k = k) + 
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

# mixed model (kept here for check only - not shown in main text. results 
# very similar to the gam but with several issues : convergence issues with 
# "cluster" effects on splines, also convergence issues on the GLM with unscaled
# variables (see below), complicated computation of odds ratios with custom
# increments. Also the readers might want to take into account the following 
# unreviewed document for a cautionary note about the interpretation of mixed
# non-gaussian GLM: 
# https://www.unige.ch/cisa/files/2917/0170/1919/CISA_BM_statsupport_20231129_GLMMcaution_verB.pdf)

# abrbin.glob.gamm <-
  # gam(
  #  conso.bin ~ s(elev.buff, k = k) + 
  #    s(northness.buff, k = k) +
  #    s(Lvrm.buff, k = k) + 
  #    s(Ldist.linear, k = k) + 
  #    s(LTirs, k = k) + 
  #    s(Appetence_mean, k = k) + 
  #    s(vis.buff,k = k)+
  #    s(Annee, k = k) + 
  #    s(Lpresence, k = k)+
  #    Bois + 
  #    Massif ,
  #  family = binomial,
  #  data = df.vrst
  # )

## Generalized Additive Model with site-specific responses ---------------------

abrbin.glob.site.gam <-
  gam(
    conso.bin ~ s(elev.buff, k = k,by = Massif) + 
      s(northness.buff, k = k,by = Massif) +
      s(Lvrm.buff, k = k,by = Massif) + 
      s(Ldist.linear, k = k,by = Massif) + 
      s(LTirs, k = k,by = Massif) + 
      s(Appetence_mean, k = k,by = Massif) + 
      s(vis.buff,k = k,by = Massif)+
      s(Annee, k = k,by = Massif) + 
      s(Lpresence, k = k,by = Massif)+
      Bois + 
      Massif ,
    family = binomial,
    data = df.vrst
  )

# check GAM residuals 
gam.check(abrbin.glob.site.gam)
summary(abrbin.glob.site.gam)

## Generalized Linear Model on scaled variables --------------------------------

abrbin.glob.lm <-
  glm(
    conso.bin ~ scale(elev.buff) + 
      scale(northness.buff) +
      scale(Lvrm.buff) + 
      scale(Ldist.linear) + 
      scale(LTirs) + 
      scale(Appetence_mean) + 
      scale(vis.buff) + 
      scale(Annee) +
      scale(Lpresence)+
      Bois + 
      Massif,
    family = binomial,
    data = df.vrst
  )

## GLM on non-scaled variables -------------------------------------------------

abrbin.nosc.glob.lm <-
  glm(
    conso.bin ~ elev.buff + 
      northness.buff +
      Lvrm.buff + 
      Ldist.linear + 
      LTirs + 
      Appetence_mean + 
      vis.buff + 
      Annee +
      Lpresence +
      Bois + 
      Massif,
    family = binomial,
    data = df.vrst
  )

# store coefficients

sc.coefs.lm <- names(fixef(abrbin.glob.lm)[2:10])

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

# note : to use mapview, linux users will have to define the web browser ahead : options(browser = 'firefox')

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
    Lvrm.buff = -7.6,
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
    "lvrm",
    "ldist.lin",
    "ltirs",
    "app_mean",
    "lvisi",
    "annee",
    "sr",
    levels(df.vrst$Bois)[2],
    levels(df.vrst$Bois)[3],
    levels(df.vrst$Massif)[2],
    levels(df.vrst$Massif)[3]
  )

# Odds ratio table for export

write.table(ci.all,"outputs/table_odds_ratios_community.txt",sep=";")

## Figure 3 : Partial residuals plots of GAM models ----------------------------

# one section per variable, then compilation of the panel. 
# includes plots for all clusters (mixed model, main text) and per clusters (SM)

# elevation 

p.elev <-
  draw(
    abrbin.glob.gam,
    residuals = T,
    select = 1,
    resid_col = "gray30",
    ci_col = "steelblue",
    smooth_col = "darkblue",
    caption = ""
  ) +
  ylim(-15,15)+
  labs(x = "Elevation (m)", "Partial effect", title = "") +
  theme_classic()

# elevation (site - specific)

p.elev.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(1, 2, 3),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(-15,15)+
  theme_classic() +
  labs(x = "Elevation (m)", "Partial effect", title = "") +
  theme_classic()

# northness 

p.north <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 2,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  ylim(-15,15)+
  labs(x = "Northness (relative index)", "Partial effect", title = "") +
  theme_classic()

# northness (site - specific)

p.north.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(4, 5, 6),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  theme_classic() +
  ylim(-15,15)+
  labs(x = "Northness (relative index)", "Partial effect", title = "") +
  theme_classic()

# rugosity

p.rugo <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 3,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  ylim(-15,15)+
  labs(x = "Rugosity (relative index, log)", "Partial effect", title = "") +
  theme_classic()

# rugosity (site specific)

p.rugo.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(7, 8, 9),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  theme_classic() +
  ylim(-15,15)+
  labs(x = "Rugosity (relative index, log)", "Partial effect", title = "") +
  theme_classic()

# distance to linear elements 

p.dist <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 4,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  labs(x = "Distance to linear element (m, log)", "Partial effect", title =
         "") +
  ylim(-15,15)+
  theme_classic()

# distance to linear elements (site specific)

p.dist.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(10, 11, 12),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  theme_classic() +
  ylim(-15,15)+
  labs(x = "Distance to linear element (m, log)", "Partial effect", title = "") +
  theme_classic()

# hunting pressure

p.hunt <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 5,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  labs(x = "Hunting shots (log)", "Partial effect", title =
         "") +
  ylim(-15,15)+
  theme_classic()

# hunting pressure (site specific)

p.hunt.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(13, 14, 15),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  theme_classic() +
  labs(x = "Hunting shots (log)", "Partial effect", title =
         "") +
  ylim(-15,15)+
  theme_classic()

# appetency

p.app <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 6,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  ylim(-15,15)+
  labs(x = "Mean appetency, relative index", "Partial effect", title = "") +
  theme_classic()

# appetency (site specific)

p.app.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(16, 17, 18),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(-15,15)+
  theme_classic() +
  labs(x = "Mean appetency", "Partial effect", title = "") +
  theme_classic()

# visibility 

p.viz <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 7,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  ylim(-15,15)+
  labs(x = "Visibility (nb pixels)", "Partial effect", title = "") +
  theme_classic()

# visibility (site specific)

p.viz.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(19, 20, 21),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(-15,15)+
  theme_classic() +
  labs(x = "Visibility (nb pixels)", "Partial effect", title = "")  +
  theme_classic()

# year effect 

p.year <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 8,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption =""
) +
  ylim(-15,15)+
  labs(x = "Years", "Partial effect", title = "") +
  theme_classic()

# year (site - specific)

p.year.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(22, 23, 24),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(-15,15)+
  theme_classic() +
  labs(x = "Years", "Partial effect", title = "")  +
  theme_classic()

# species richness

p.sr <- draw(
  abrbin.glob.gam,
  residuals = T,
  select = 9,
  resid_col = "gray30",
  ci_col = "steelblue",
  smooth_col = "darkblue",
  caption = ""
) +
  ylim(-15,15)+
  labs(x = "Species richness (log)", "Partial effect", title =
         "") +
  theme_classic()

# species richness (site-specific)

p.sr.site <- abrbin.glob.site.gam |>
  draw(grouped_by = TRUE,
       select = c(25, 26, 27),
       caption = "") +
  scale_color_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  scale_fill_manual(
    values = c('#46327e', '#1fa187', 'goldenrod'),
    name = "cluster",
    labels = c("Semnoz", "Cimeteret", "Hautes Bauges")
  ) +
  ylim(-15,15)+
  theme_classic() +
  labs(x = "Species richness (log)", y = "Partial effect", title = "")

# Categorical variables

labs.type <- c("young", "medium", "old")
p.stype <- parametric_effects(abrbin.glob.gam, terms ='Bois') |> 
draw(caption = "") +
theme_classic() +
  labs(x = "Stand type",y = "Partial effect", title = "",caption = "") +
  ylim(-0.5,0.5)+
 scale_x_discrete(labels= labs.type)  

labs.cluster <- c("Semnoz","Cimeteret","Hautes Bauges")
p.cluster <- parametric_effects(abrbin.glob.gam, terms ='Massif') |> 
  draw(caption = "") +
  theme_classic() +
  scale_x_discrete(labels= labs.cluster)  +
  ylim(-0.5,0.5)+
  labs(x = "Cluster",y = "Partial effect", title = "",caption = "") 

cowplot::plot_grid(
  p.sr,
  p.app,
  
  p.elev,
  p.north,
  p.rugo,
  
  p.dist,
  p.hunt,
  p.viz,
  
  p.year,
  p.stype,
  p.cluster,
  
  nrow = 4,
  ncol = 3,
  labels = "auto"
)

ggsave("outputs/gaudry_et_al_fig3.svg",
       height = 12,
       width = 9)

ggsave("outputs/gaudry_et_al_fig3.png",
       height = 12,
       width = 9)

## MULTISPECIES COMPARATIVE ANALYSIS ===========================================

## Table 2 : descriptive data on plant species ---------------------------------

abr.sp.mean <- aggregate(prop.abr[,c("Presence","Consommation")],by = list(prop.abr$Nom_Latin), FUN = "mean")
abr.sp.sd <- aggregate(prop.abr[,c("Presence","Consommation")],by = list(prop.abr$Nom_Latin), FUN = "sd")
abr.sp.rg <- aggregate(prop.abr[,c("Presence","Consommation")],by = list(prop.abr$Nom_Latin), FUN = "range")

tab.2 <- cbind(abr.sp.mean,abr.sp.sd[,-1],abr.sp.rg$Presence,abr.sp.rg$Consommation)
colnames(tab.2) <- c("species","mean_occ","mean_browsed","sd_occ","sd_browsed","min_occ","max_occ","min_browsed","max_browsed")
tab.2[,-1] <- round(tab.2[,-1],2)
write.csv(tab.2,"outputs/Table2.csv")

## yearly trends in browsing per species ---------------------------------------

#on plots with at least one species browsed

df.sp.trends <- data.frame()

for(i in 1:length(sp.keep)){
  sp.test <- sp.keep[i]
  tp.trend <- subset(prop.abr, Nom_Latin == sp.test)
  off.trend <- log(tp.trend$Presence)
  glm.trend <- glm(Consommation ~ Annee, offset = off.trend,data = tp.trend, family = poisson)
  
  coefs.trend <- summary(glm.trend)$coefficients
  
  pval.over <- testDispersion(glm.trend)$p.value
  
  df.sp.trends <- rbind(df.sp.trends,c(coefs.trend[2,],pval.over))

}
colnames(df.sp.trends) <- c("trend","se.trend","z.trend","pval.trend","overdisp.trend")
rownames(df.sp.trends) <- sp.keep

# check Abies alba (only significant linear trend among all species)

tp.trend.aa <- subset(prop.abr, Nom_Latin == "Abies alba")
off.trend.aa <- log(tp.trend.aa$Presence)
glm.trend.aa <- glm(Consommation ~ Annee, offset = off.trend.aa,data = tp.trend.aa, family = poisson)
summary(glm.trend.aa)

## prepare dataset for species-level models ------------------------------------

covariates <- df.vrst[, !names(df.vrst) %in% "conso.bin"]

abr.data4 <-
  merge(abr.data3,
        covariates,
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

## multi-species mixed GAM (not used in the article) ---------------------------

#abrbin.sp.gamm <-
"gamm(
    Consommation ~ s(elev.buff, k = k, by = Nom_Latin) + 
      s(northness.buff, k = k, by = Nom_Latin) +
      s(Lvrm.buff, k = k, by = Nom_Latin) + 
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
             s(Lvrm.buff, k = k) + 
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
             scale(Lvrm.buff) + 
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
             Lvrm.buff + 
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

## comparison of species responses to each variable  ---------------------------

# compare the curves(uses gratia) - for exploration only, not used in article

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
      "s(Lvrm.buff)",
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
# the order of levels in the "model" variable is forced, otherwise unnest reorders randomly... 

comp.sp.all.df <- unnest(comp.sp.all,cols="data")
comp.sp.all.df$model <- factor(comp.sp.all.df$.model , levels = unique(comp.sp.all.df$.model)) 

p.gam.ele <- ggplot(comp.sp.all.df) +
  aes(x = elev.buff, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line(linewidth = 2) +
  theme_classic() +
  labs(x = "Elevation (m)", y = "Scaled marginal effect", title =
         "") 

p.gam.north <- ggplot(comp.sp.all.df) +
  aes(x = northness.buff, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,3) + 
  labs(x = "Northness (relative unit)", y = "Scaled marginal effect", title =
         "")
p.gam.vrm <- ggplot(comp.sp.all.df) +
  aes(x = Lvrm.buff, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,2) +
  labs(x = "Rugosity (log(relative unit))", y = "Scaled marginal effect", title =
         "")

p.gam.lin <- ggplot(comp.sp.all.df) +
  aes(x = Ldist.linear, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "Distance to near.estimate linear element (log(m))", y = "Scaled marginal effect", title =
         "")

p.gam.tirs <- ggplot(comp.sp.all.df) +
  aes(x = LTirs, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "number of hunting shots (log)", y = "Scaled marginal effect", title =
         "")

p.gam.app <- ggplot(comp.sp.all.df) +
  aes(x = Appetence_mean, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "mean Appetence", y = "Scaled marginal effect", title =
         "") 

p.gam.vis <- ggplot(comp.sp.all.df) +
  aes(x = vis.buff, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1,1)+
  labs(x = "Visbility (Unit???)", y = "Scaled marginal effect", title =
         "") 

p.gam.an <- ggplot(comp.sp.all.df) +
  aes(x = Annee, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-2, 2) +
  labs(x = "Years", y = "Scaled marginal effect", title =
         "")

p.gam.sr <- ggplot(comp.sp.all.df) +
  aes(x = Lpresence, y = .estimate, col = .model) +
  scale_color_viridis_d(labels = levels(abr.data5$Nom_Latin)) +
  geom_line() +
  theme_classic() +
  ylim(-1, 3) +
  labs(x = "Species richness", y = "Scaled marginal effect", title =
         "")
# figure

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

# ggsave("outputs/species_level_plots.png", width = 9, height = 9)

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
  
  # compile data
  
  ci.sp <- as.data.frame(elev.sp.or[, c(2, 3, 4)])
  colnames(ci.sp) <- c("odds", "lower", "upper")
  ci.sp$species <-
    sub("_", " ", substring(glm.sp.list[i], 11, nchar(glm.sp.list[i])))
  ci.sp$variable <-
    c(
      "elev",
      "northness",
      "lvrm",
      "ldist.lin",
      "ltirs",
      "app_mean",
      "visi",
      "annee",
      "sr",
      levels(abr.data5$Bois)[2],
      levels(abr.data5$Bois)[3],
      levels(abr.data5$Massif)[2],
      levels(abr.data5$Massif)[3]
    )
  or.glm <- rbind(or.glm, ci.sp)
}

# concatenate with the odds ratios of the community model

ci.sp.all <- rbind(or.glm,ci.all)
write.csv2(ci.sp.all,"outputs/odds_ratios_glm.csv")

# plots odds ratios as a forest plot

varinames <- data.frame(names = c("elev","northness","lvrm","ldist.lin","ltirs","app_mean","annee","sr"),
                        varilabs = c("Elevation (m)","Northness (relative unit","Rugosity (log(relative unit))","Distance to nearest linear element (log(m))","Number of hunting shots (log)","Mean attractivity (scores)","Years","Species richness")
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

# plot of "stand" effect (contrasts, reference = "PB" (young stands))

or.bois.sp <- subset(ci.sp.all, variable%in%c(levels(abr.data5$Bois)[2],levels(abr.data5$Bois)[3])& species!="community")
or.bois.sp$species <- reorder(or.bois.sp$species,or.bois.sp$odds)
or.bois.sp$n_species <- as.numeric(or.bois.sp$species)

or.bois.comm <- subset(ci.sp.all, variable%in%c(levels(abr.data5$Bois)[2],levels(abr.data5$Bois)[3])& species=="community")
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

# plot of the cluster effect (contrasts, reference = "Semnoz")

or.zone.sp <- subset(ci.sp.all, variable%in%c(levels(abr.data5$Massif)[2],levels(abr.data5$Massif)[3])& species!="community")
or.zone.sp$species <- reorder(or.zone.sp$species,or.zone.sp$odds)
or.zone.sp$n_species <- as.numeric(or.zone.sp$species)

or.zone.comm <- subset(ci.sp.all, variable%in%c(levels(abr.data5$Massif)[2],levels(abr.data5$Massif)[3])& species=="community")
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

## Figure 4a: forest plot of scaled effects - community level ------------------

labs.variables <- c("log(Species richness)",
                    "Appetency",
                    "Years",
                    "Elevation",
                    "log(Hunting shots)",
                    "Distance",
                    "log(Rugosity)",
                    "Visibility",
                    "Northness")

# plot - level model (the species-level plots are in a next section)

sc.coefs.lm <- names(coef(abrbin.glob.lm)[2:10])

fp.comm <- plot_model(
  abrbin.glob.lm,
  terms = sc.coefs.lm,
  sort.est = T,
  transform = NULL,
  axis.lim = c(-2, 2),
  colors = c("steelblue","darkred"),
  dot.size = 2, 
  line.size = 0.5,
  vline.color = "gray80"
) +
  theme_classic() +
  labs(x = "", y = "Estimates - 95%CI", title = paste(letters[1],"Plot", sep = " - "))+
  scale_x_discrete(
    labels =
      c(
        "Northness",
        "Visibility",
        "log(Rugosity)",
        "Distance",
        "log(Hunting shots)",
        "Elevation",
        "Years",
        "Appetency",
        "log(Species richness)"
      )
  )

fp.comm



# ggsave("outputs/fit_scaled_global.png",height=8,width=4)

# same thing pour non scaled coefficients 

nosc.coefs.lm <- names(coef(abrbin.nosc.glob.lm)[2:10])

fp.comm.nosc <- plot_model(
  abrbin.nosc.glob.lm,
  terms = nosc.coefs.lm,
  colors = c("steelblue","darkred"),
  sort.est = T,
  transform = NULL, 
  show.values = T,
  value.offset = 0.3,
  value.size = 2,
  dot.size = 1, 
  line.size = 0.5,
  vline.color = "gray80"
) +
  theme_classic() +
  labs(x = "", y = "Estimates - 95%CI", title = paste(letters[1],"Plot", sep = " - "))+
  scale_x_discrete(labels=rev(      c(
    "Northness",
    "Visibility",
    "log(Rugosity)",
    "Distance",
    "log(Hunting shots)",
    "Elevation",
    "Years",
    "Appetency",
    "log(Species richness)"
  ))+
                     theme( axis.title.x = element_text(size = 12))
)

fp.comm.nosc
# ggsave("outputs/fit_scaled_global_nosc.png",height=8,width=4)
 

## Figure 4b: forest plot of scaled effects - species level --------------------

# has to be combined with the first part of figure 4 above
# order plots as in the plot level analysis

names(glm.sp.list) <- levels(abr.data5$Nom_Latin)
coef.plot.level <- names(rev(sort(coef(abrbin.glob.lm)[2:10])))


for (i in 1:length(glm.sp.list)){
  
glm.tp <- get(glm.sp.list[i])
sc.coefs.lm <- names(coef(glm.tp)[2:10])

pos.coefs <- match( coef.plot.level,sc.coefs.lm)

if(i %in%c(1,2,4,5,7,8,10)){
assign(paste("sp_forplot",i,sep="_"),plot_model(
  glm.tp,
  order.terms = pos.coefs,
  terms = coef.plot.level,
  sort.est = F,
  colors = c("steelblue","darkred"),
  transform = NULL,
  axis.labels = rep(" ",9),
  axis.lim = c(-2,2),
 dot.size = 1, 
 line.size = 0.5,
 vline.color = "gray80"
) +
  theme_classic() +
  labs(x = "", y = "Estimates - 95%CI", title = paste(letters[i+1],names(glm.sp.list)[i], sep = " - "))
)
} else {
  
  assign(paste("sp_forplot",i,sep="_"),plot_model(
    glm.tp,
    terms = coef.plot.level,
    sort.est = F,
    order.terms = pos.coefs,
    colors = c("steelblue","darkred"),
    transform = NULL,
    axis.lim = c(-2,2), dot.size = 1, 
    line.size = 0.3,
    vline.color = "gray80"
    ) +
    theme_classic() +
    labs(x = "", y = "Estimates - 95%CI", title = paste(letters[i+1],names(glm.sp.list)[i], sep = " - "))+
      scale_x_discrete(labels=rev(labs.variables))
  )
}

}

fp.comm+ 
sp_forplot_1+
sp_forplot_2+
sp_forplot_3+
sp_forplot_4+
sp_forplot_5+
sp_forplot_6+
sp_forplot_7+
sp_forplot_8+
sp_forplot_9+
sp_forplot_10+
  plot_layout(ncol = 3,widths = 1,heights = 1)

ggsave("outputs/species_level_forest_plots.png", width = 12, height = 12)

## Appendix S1 : data for maps of covariates -----------------------------------

# main layer 
sm1 <- as.data.frame(spat.data3)
sm1.sf <- st_as_sf(sm1)
#st_write(sm1.sf, dsn = "outputs/SM1.shp")

# Bauges Regional Park
bnrp <- st_read(dsn = "data/PNR_bauges.geojson")

# city markers
cities <- st_read(dsn = "data/gaudry_et_al_city_markers.geojson")

# France - contour
fra <- st_read(dsn = "data/France_contour.geojson")

tm_shape(fra)+
  tm_borders() 

# map boundaries
map.box <- st_bbox(bnrp)
xrange <- map.box$xmax - map.box$xmin 
yrange <- map.box$ymax - map.box$ymin 

map.box[1] <- map.box[1] - (0.1 * xrange) 
map.box[3] <- map.box[3] + (0.1 * xrange) 
map.box[2] <- map.box[2] - (0.1 * yrange) 
map.box[4] <- map.box[4] + (0.1 * yrange) 
map.box <- map.box %>%  
  st_as_sfc() 

# baseline map
base.sm1 <- tm_shape(bnrp,bbox = map.box)+
  tm_fill(col = "#1e9b8a",alpha = 0.2)+
  tm_shape(cities) + 
  tm_markers(size = 0.4, col = "black",shape = 17,text = "nom_commune",text.size = 0.7,xmod = 2)+
  #tm_symbols(size = 0.4, col = "black",shape = 17)+
  tm_compass(type = "arrow", position = c("right", "bottom"),show.labels = 0,color.dark = "black",color.light = "white",lwd=1)+ 
  tm_scale_bar(breaks = c(0, 5, 10),  position = c("right", "bottom"),lwd = 0.5,text.size = 0.7)

# loop over all variables 

labs.map <-
  c(
    "Elevation (m)",
    "Slope",
    "Northness (relative units)",
    "Rugosity (relative units)",
    "Distance to \n linear elements (m)",
    "Human frequentation",
    "Hunting shots (count)",
    "Visibility (pixel count)"
  )

names(labs.map) <- c("elevation","slope","northness","rugosity","distance","human","hunting","visibility")

for(i in 1:8){
  
var <- colnames(sm1.sf)[i+2]
map.var <- base.sm1+
  tm_shape(sm1.sf)+
  tm_dots(col = var,size = 0.2,shape = 21,palette = "viridis",title = labs.map[i],title.size = 0.7,style = "pretty")

nm.map <- paste("sm2",colnames(sm1.sf)[i+2],sep=".")

assign(nm.map,map.var)
 tmap_save(map.var, paste("outputs/S1_",names(labs.map[i]),".png",sep=""), dpi = 600)

}


## SM 2 : list of plant species ------------------------------------------------

abr.data3b <- subset(abr.data2b, conso.bin == 1)

# number of plots of occurrence per year

occ.per.year <-
  aggregate(
    abr.data3b[,c("Presence","Consommation")],
    by = list(abr.data3b$Nom_Latin, abr.data3b$Annee),
    FUN = function(x){sum(x>0)}
  )
colnames(occ.per.year)[1:2] <- c("species","year")

abr.sp.mean.all <-
  aggregate(occ.per.year[, c("Presence", "Consommation")],
            by = list(occ.per.year$species),
            FUN = "mean")

abr.sp.sd.all <-
  aggregate(occ.per.year[, c("Presence", "Consommation")],
            by = list(occ.per.year$species),
            FUN = "sd")

abr.sp.rg.all <-
  aggregate(occ.per.year[, c("Presence", "Consommation")],
            by = list(occ.per.year$species),
            FUN = "range")

tab.2.all <-
  cbind(abr.sp.mean.all,
        abr.sp.sd.all[, -1],
        abr.sp.rg.all$Presence,
        abr.sp.rg.all$Consommation)

colnames(tab.2.all) <-
  c(
    "species",
    "mean_occ",
    "mean_browsed",
    "sd_occ",
    "sd_browsed",
    "min_occ",
    "max_occ",
    "min_browsed",
    "max_browsed"
  )

tab.2.all[,-1] <- round(tab.2.all[,-1],2)

tab.sm2 <-
  merge(tab.2.all,
        app.sp,
        by.x = "species",
        by.y = "espece",
        all.x = T, all.y = F)

#write.xlsx(tab.sm2, "outputs/SM2.xlsx", sheetName = "SM2")


## SM 3 : model coefficients and odds ratios -----------------------------------

## coefficients of qualitative variables

# summarises all coefficients as a table - scaled variables

eval.glm.sp <- lapply(glm.sp.list, FUN = "get")
eval.glm.all <- c(list(abrbin.glob.lm), eval.glm.sp)

tab.scaled <-
  modelsummary(
    eval.glm.all,
    estimate = "{estimate} [{conf.low}, {conf.high}]",
    statistic = NULL,
    shape = model + statistic ~ term,
    fmit = fmt_significant(3),
    output = "outputs/glm-scaled-variables.xlsx"
  )

# summarises all coefficients as a table - non scaled variables

eval.glm.nosc.sp <- lapply(glm.nosc.sp.list, FUN = "get")
eval.glm.nosc.all <- c(list(abrbin.nosc.glob.lm), eval.glm.nosc.sp)

tab.raw <-
  modelsummary(
    eval.glm.nosc.all,
    estimate = "{estimate} [{conf.low}, {conf.high}]",
    statistic = NULL,
    shape = model + statistic ~ term,
    fmit = fmt_significant(3),
    output = "outputs/glm-raw-variables.xlsx"
  )

# odds ratios


# creates the main table
write.xlsx(tab.scaled, "outputs/SM3.xlsx", sheetName = "scaled_coef")
write.xlsx(tab.raw, "outputs/SM3.xlsx", sheetName = "raw_coef")


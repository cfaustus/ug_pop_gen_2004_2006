##########
# Code for population genetic diversity metrics
# associated with Faust et al. 2019
# Parasites & Vectors:
# run with R version 3.6.1 (2019-07-05)

##########
# libraries
library(adegenet)
library(poppr) 
library(pegas)
library(mmod)
library(ape)
library(hierfstat)
library(ggplot2)

##########
# data import
schisto_full = read.genepop("data/microsat.gen")
strata = read.csv("data/strata_per_pop.csv",header=TRUE) #change to metadata in enlighten
head(strata)
schistodf = genind2df(schisto_full,sep="/") #convert to a dataframe
head(schistodf)
pop = seppop(schisto_full)
schisto_full = df2genind(schistodf[, -1], sep = "/", pop = schistodf$pop, strata = strata)
schisto_full

schisto_trun = missingno(schisto_full, type="geno",cutoff=0) # Remove missing data using poppr

########
# number of alleles per locus
locus_table(schisto_trun)

########
# Population summary metrics (poppr)
class(child_wk)
# first define population, then output is a popprtable (also data frame)
# Columns are as follows:
# N:	Number of miracidia observed.
# MLG: Number of multilocus genotypes (MLG) observed.
# eMLG: Number of expected MLG, based on rarefaction
# SE: Standard error based on eMLG.
# H: MLG diversity Shannon-Wiener Index [Shannon 2001]
# G: MLG diversity Stoddart and Taylor’s Index [Stoddart & Taylor 1988]
# lambda:	Simpson’s Index [Simpson 1949]
# E.5:	Evenness/ E5 [Pielou 1975; Ludwig & Reynolds 1988; Grünwald et al. 2003]
# Hexp:	Nei’s unbiased gene diversity [Nei 1978]
# Ia: Index of association, IA [Brown, Feldman & Nevo 1980; Smith et al. 1993]
# rbarD: Standardized index of association

setPop(schisto_trun) = ~ ChildID/Week
child_wk = poppr(schisto_trun)
setPop(schisto_trun) = ~ ChildID #aka infrapopulation
infrapop = poppr(schisto_trun)
setPop(schisto_trun) = ~ School
school = poppr(schisto_trun)
setPop(schisto_trun) = ~ School/Week
school_wk = poppr(schisto_trun)

#######
# allelic richness, observed heterozygousity and exp heterozygousity
# per child per timepoint
setPop(schisto_trun) = ~ ChildID/Week
dat_hfs = genind2hierfstat(schisto_trun) # converting to hierfstat format
populations = unique(dat_hfs$pop)
ar = as.data.frame(allelic.richness(dat_hfs)) # outputs allelic richness; rows = alleles; col = AR
alleles = rownames(ar)
ar = t(ar)
colnames(ar) = alleles
ar = ar[-1,]
ar_sum = as.data.frame(matrix(NA, nrow = nrow(ar)))
mean = as.vector(rowMeans(ar[,1:7]))
ar_sum$ar_mean = mean
ar_sum$pop = populations
dat_hfs_pop = seppop(schisto_trun) 
summary_by_pop = lapply(dat_hfs_pop, summary) 
summary_df = as.data.frame(matrix(NA, nrow = length(summary_by_pop), ncol = 3)) 
for (i in 1:length(summary_by_pop)){ 
  summary_df[i,1] = names(summary_by_pop[[i]]$pop)
  summary_df[i,2] = summary_by_pop[[i]]$n.by.pop
  summary_df[i,3] = mean(summary_by_pop[[i]]$Hobs)
  summary_df[i,4] = mean(summary_by_pop[[i]]$Hexp)
} 
colnames(summary_df) = c('pop', 'n_mir', 'Hobs', 'Hexp')
div_stat = merge(summary_df, ar_sum) 
div_stat = div_stat[ , -which(names(div_stat) %in% c("V1"))]
#write.csv(div_stat,"output/diversitystats.csv", row.names =F)

#######
# AMOVA on subset 
setPop(schisto_trun) = ~ Week #set population by week
wkind = unique(strata$Week) # number of unique weeks
amova_df = as.data.frame(matrix(data=NA, 
                                nrow=length(wkind)*3,ncol=8)) # empty dataframe
colnames(amova_df) = c('week', 'level', 'Df', 'Sum Sq', 
                       'Mean Sq','Sigma', '%', 'p_value')
for (i in 1:length(wkind)){
  amv_pop = popsub(schisto_trun, paste(wkind[i]))
  schisto_amova = poppr.amova(amv_pop, ~School/ChildID, within = FALSE)
  set.seed(1999)
  schistosignif   = randtest(schisto_amova, nrepet = 1000)
  plot(schistosignif)
  amova_df[((i-1)*4+1):((i-1)*4+4),'week'] = wkind[i]
  amova_df[((i-1)*4+1):((i-1)*4+4),'level'] = row.names(schisto_amova$results)
  amova_df[((i-1)*4+1):((i-1)*4+4),'Df'] = schisto_amova$results[['Df']]
  amova_df[((i-1)*4+1):((i-1)*4+4),'Sum Sq'] = schisto_amova$results[['Sum Sq']]
  amova_df[((i-1)*4+1):((i-1)*4+4),'Mean Sq'] = schisto_amova$results[['Mean Sq']]
  amova_df[((i-1)*4+1):((i-1)*4+4),'Sigma'] = schisto_amova$componentsofcovariance[['Sigma']]
  amova_df[((i-1)*4+1):((i-1)*4+4),'%'] = schisto_amova$componentsofcovariance[['%']]
  amova_df[((i-1)*4+1):((i-1)*4+3),'p_value'] = schistosignif$pvalue
}
#write.csv(amova_df, 'output/amova_trun_data.csv', row.names = F)
# uncomment to write output file

### cumulative population amova
schisto_trun_amova = poppr.amova(schisto_trun, ~School/ChildID, within = FALSE)
schisto_signif   = randtest(schisto_trun_amova, nrepet = 1000) #takes a long time to run

#########
# FST
Gst_Hedrick(schisto_trun) # Hendrick’s standardized GST
setPop(schisto_trun) = ~ Week/ChildID
fstat(schisto_trun, fstonly = TRUE)

setPop(schisto_trun) = ~ Week/School
wksch_fst = pairwise.fst(schisto_trun, res.type="matrix")
wksch_tree = nj(wksch_fst)
plot(wksch_tree, type="unr", font=2)
annot = round(wksch_tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
table.paint(wksch_fst, col.labels=1:16)
wksch_fst_plot = wksch_fst
diag(wksch_fst_plot) = NA
sch_time = row.names(wksch_fst)
wksch_fst_plot = wksch_fst_plot[order(factor(sch_time, levels = c("0_BUG","4_BUG","26_BUG","27_BUG","52_BUG",
                                      "53_BUG","56_BUG","104_BUG","105_BUG","108_BUG",
                                      "0_BW","1_BW","4_BW","26_BW","27_BW",     
                                      "52_BW","53_BW", "104_BW", "105_BW", "108_BW",
                                      "0_MUS", "1_MUS", "4_MUS", "26_MUS", "27_MUS",
                                      "52_MUS","53_MUS", "56_MUS","104_MUS",
                                      "105_MUS","108_MUS"))), 
    order(factor(sch_time, levels = c("0_BUG","4_BUG","26_BUG","27_BUG","52_BUG",
                                      "53_BUG","56_BUG","104_BUG","105_BUG","108_BUG",
                                      "0_BW","1_BW","4_BW","26_BW","27_BW",     
                                      "52_BW","53_BW", "104_BW", "105_BW", "108_BW",
                                      "0_MUS", "1_MUS", "4_MUS", "26_MUS", "27_MUS",
                                      "52_MUS","53_MUS", "56_MUS","104_MUS",
                                      "105_MUS","108_MUS")))]
boxplot(wksch_fst_plot,  xlab="Population", 
        ylab="Fst", las = 2, lab.cex = 0.2)
#write.csv(wksch_fst_plot, "output/fst_sch_time.csv", row.names = TRUE)

########
# Trees of repeatedly sampled infrapopulations:
# BUG30619; MUS20615; MUS20625; MUS20626; MUS30622; MUS807
setPop(schisto_trun) = ~ ChildID

MUS807 = popsub(schisto_trun, sublist = "MUS807",drop = TRUE) #can change to others
MUS807 %>%
  genind2genpop(pop = ~Week/School) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)


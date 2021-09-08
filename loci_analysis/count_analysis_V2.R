#count_analysis_v2 010621

library(tidyverse)
library(data.table)
library(viridis)
library(grid)
library(UpSetR)
library(plot3D)

## Read in data from various sources
retroseq.matrix <- read_csv("0602_fil_retroseq_locifl7q60_pad100_countfl5q30_callmatrix.csv", col_names = T) %>%
  separate(col = X1, into = c("Chromosome", "Coordinate"), sep = "-") %>% 
  mutate(method = "retroseq") %>% 
  relocate(method)

cnv.matrix <- read_table2("1004_cnvnator_indmatrix.txt", col_names = T) %>% 
  separate(col = locus, into = c("method", "Chromosome", "Coordinate"), sep = "_") %>% 
  mutate(Coordinate = gsub("-","_",Coordinate))

delly.matrix <- read_table2("1004_delly_indmatrix.txt", col_names = T) %>% 
  separate(col = locus, into = c("method", "Chromosome", "Start", "End"), sep = "_") %>%
  unite("Coordinate", Start:End, sep = "_") %>% 
  mutate(method = "delly") %>% 
  rename_with(~ stringr::str_replace_all(.,"_", "\\."))


## Check for matching column names
all_equal(retroseq.matrix, cnv.matrix) # Different number of columns
setdiff(colnames(retroseq.matrix), colnames(cnv.matrix)) #fortisxfuliginosa.Daphne.14612
setdiff(colnames(retroseq.matrix), colnames(delly.matrix)) #fortisxfuliginosa.Daphne.14612
# One sample is unique to retroseq. Will DROP it

## Melt into long format
##Combine into single tibble with loci present based on count column. cnv and delly inverted

retroseq.df <- retroseq.matrix %>%
  gather("sampleID","count",4:296) %>%
  filter(sampleID != "fortisxfuliginosa.Daphne.14612") %>% 
  separate(sampleID, into = c("species","island","idNumber")) %>% 
  mutate(locusPresent = if_else((count < 1) , 0, 1))

cnv.df <- cnv.matrix %>%
  gather("sampleID","count",4:295) %>%
  separate(sampleID, into = c("species","island","idNumber")) %>% 
  mutate(locusPresent = if_else((count == 1) , 0, 1)) %>% 
  mutate(method = "cnv")

delly.df <- delly.matrix %>%
  gather("sampleID","count",4:295) %>%
  separate(sampleID, into = c("species","island","idNumber")) %>% 
  mutate(locusPresent = if_else((count >= 2) , 0, 1)) %>% 
  mutate(locusPresent = replace_na(locusPresent, 1))

erv_counts.df <- bind_rows(retroseq.df, cnv.df) %>%
  bind_rows(delly.df) %>% 
  mutate(assembly = if_else(grepl("chr",Chromosome),"chromosome","scaffold")) %>%
  mutate(Coordinate = gsub("_","-",Coordinate)) %>% 
  mutate(locusID = paste(Chromosome, Coordinate, sep = "_")) %>% 
  # filter(assembly == "chromosome") %>%
  select(-count)

#Check for NAs
retroseq.df %>% 
  filter_all(any_vars(is.na(.)))
cnv.df %>% 
  filter_all(any_vars(is.na(.)))
delly.df %>% 
  filter_all(any_vars(is.na(.))) # Many counts are NA, Mette advised keep as present loci
erv_counts.df %>% 
  filter_all(any_vars(is.na(.))) # No NAs

# save(erv_counts.df, file = "erv_counts.df.210618.Rdata")

load("n.all.Rdata")
# load("erv_counts.df.210618.Rdata")

## Save method used for each locus for later annotation
Loci_Method.df <- erv_counts.df %>% 
  select(locusID, method) %>% 
  unique()

spplist <- c("noctis" , "bicolor", "olivacea","fusca","crassirostris", 
             "inornata", "parvulus",
             "pauper","heliobates",
             "pallidus","psittacula",
             "difficilis", "fuliginosa" ,"acutirostris", "BB", 
             "fortis","fortisxfuliginosa","forxfuligxscandens", "forxsca", "conirostris",
             "magnirostris","scandens", "septentrionalis","propinqua")
spp.list <- list()
spp.list$allspp <- spplist
spp.list$complex <- c("parvulus","fuliginosa","fortis","scandens","fortisxfuliginosa","forxfuligxscandens","forxsca","BB")
spp.list$other <- subset(spp.list$allspp,!(spp.list$allspp %in% spp.list$complex))
spp.list$outgroups <- c("noctis","bicolor")

## Load identified ervs and loci
loci_erv_retroseq_summary <- read_csv("0602_fil_retroseq_locifl7q60_pad100_countfl5q30_ervsummary.csv")
loci_erv_delly_summary <- read_csv("1023_updated_dellysummary.csv")
loci_erv_cnv_summary <- read_csv("1025_cnv300str.blatmatch.summary.csv")

## Summarize loci
locusID_ervAssigned_RS <- loci_erv_retroseq_summary %>%
  rename_with(~ c("locus", "erv.call", "erv.assigned.freq",
                  "n.erv.called", "end", "start",
                  "all.hits.ervs","all.hits.counts","locus.width","width")) %>% 
  separate(locus, into = c("chr", "x", "y")) %>% 
  mutate(locusID = paste(chr,"_",x,"-",y, sep = "")) %>%
  mutate(erv.assigned = str_extract(erv.call, paste(n.all, collapse = "|"))) %>% 
  select(locusID, erv.assigned) %>%
  mutate(method = "retroseq")

locusID_ervAssigned_delly <- loci_erv_delly_summary %>%
  mutate(locusID = paste(chr,"_",start,"-",end, sep = "")) %>%
  mutate(erv.assigned = str_extract(erv.call, paste(n.all, collapse = "|"))) %>% 
  select(locusID, erv.assigned) %>% 
  mutate(method = "delly")

locusID_ervAssigned_cnv <- loci_erv_cnv_summary %>%
  mutate(locusID = paste(chr,"_",start,"-",end, sep = "")) %>%
  mutate(erv.assigned = str_extract(erv.call, paste(n.all, collapse = "|"))) %>% 
  select(locusID, erv.assigned) %>% 
  mutate(method = "cnv")

locusID_ervAssigned_all <- locusID_ervAssigned_RS %>% 
  bind_rows(locusID_ervAssigned_delly) %>% 
  bind_rows(locusID_ervAssigned_cnv) #21,258 loci

locusID_ervAssigned_all %>% filter_all(any_vars(is.na(.))) # No NAs

## Count erv hits to get an idea about which ervs might be attracting hits
finchERV_score_summary <- read.table("FinchERV_score_summary.tsv", 
                                     sep = "\t", 
                                     header = T, 
                                     stringsAsFactors = F) %>% 
  rename(erv.assigned = erv)

loci_count.erv_assigned.df <- erv_counts.df %>% 
  left_join(locusID_ervAssigned_all)

loci_count.erv_assigned.df %>% filter_all(any_vars(is.na(.))) # No NAs

erv_hits.summary <- loci_count.erv_assigned.df %>% 
  group_by(erv.assigned, method) %>% 
  filter(locusPresent == 1) %>% 
  tally() %>% 
  left_join(finchERV_score_summary) %>% 
  spread(method, n) %>% 
  replace_na(list(score = 2000, subgenes = "reference", delly = 0, cnv = 0, retroseq = 0)) %>%
  mutate(n = delly + cnv + retroseq)

top10ervs <- erv_hits.summary %>% 
  # filter(n >= 9400) %>% 
  filter(n >= 14400) %>%
  arrange(-n)

top5ervs <- erv_hits.summary %>% 
  filter(n >= 40000) %>% 
  arrange(-n)

## Get tips from phylo

ervTree <- ape::read.nexus("camPar1_RetroAligner_AAguided_trim10gap45.FastTree2gtr.converted.tre")

## Plot Hits by ERV in phylogenetic order

ggplot(data = erv_hits.summary %>% 
         mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
         droplevels()) +
  geom_point(aes(x = erv.assigned,
                 y = n),
             size = 1.5,
             shape = 18) +
  scale_x_discrete() +
  # scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))


## Summary with normalized count by sample size
sampling_counts.species <- erv_counts.df %>% 
  # select(species, island) %>% #Old incorrect calc
  select(species, island, idNumber) %>%
  distinct() %>% 
  group_by(species) %>% 
  summarise(birds_sampled = n())

loci_count.erv_assigned.top10.spp.summary <- loci_count.erv_assigned.df %>% 
  filter(erv.assigned %in% c(top10ervs$erv.assigned)) %>% 
  filter(locusPresent == 1) %>% 
  group_by(species, erv.assigned) %>% 
  summarise(group_counts = sum(locusPresent)) %>% 
  left_join(sampling_counts.species, by = c("species")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(ervTree$tip.label[ervTree$tip.label %in% top10ervs$erv.assigned])))

sampling_counts.species.island <- erv_counts.df %>% 
  select(species, island, idNumber) %>% 
  distinct() %>% 
  group_by(species, island) %>% 
  summarise(birds_sampled = n())

loci_count.erv_assigned.top10.spp.island.summary <- loci_count.erv_assigned.df %>% 
  filter(erv.assigned %in% c(top10ervs$erv.assigned)) %>% 
  filter(locusPresent == 1) %>% 
  group_by(species, erv.assigned, island) %>% 
  summarise(group_counts = sum(locusPresent)) %>% 
  left_join(sampling_counts.species.island, by = c("species", "island")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(ervTree$tip.label[ervTree$tip.label %in% top10ervs$erv.assigned])))

loci_count.erv_assigned.allERVs.spp.summary <- loci_count.erv_assigned.df %>% 
  filter(locusPresent == 1) %>% 
  group_by(species, erv.assigned) %>% 
  summarise(group_counts = sum(locusPresent)) %>% 
  left_join(sampling_counts.species, by = c("species")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled)

## Histogram
p <- erv_counts.df %>% 
  group_by(locusID) %>%
  summarise(locus_frequency = sum(locusPresent)/292) %>%
  left_join(locusID_ervAssigned_all) %>%
  ungroup() %>% 
  filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = top10ervs$erv.assigned)) %>% 
  ggplot(aes(x = locus_frequency)) +
  geom_histogram(bins = 50, fill = "orange", color = "black") +
  facet_wrap(~erv.assigned, scales = "free_y", ncol = 2, nrow = 5) +
  ggtitle("Frequency of top 10 erv loci in all species")

pdf(file = "plots/top10_ERV_freq_hist_210619.pdf", height = 8, width = 8)
p
dev.off()

## Heatmap

p <- ggplot(loci_count.erv_assigned.top10.spp.summary,
            aes(fill = (ervCounts_normalized), x = species, y = erv.assigned)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  # labs(fill = "ln(ERV loci)") +
  labs(fill = "ERV loci") +
  ggtitle("Top 10 ERV loci by species")

pdf(file = "plots/top10_ERV_heatmap_210630.pdf", height = 4.4, width = 3.75)
# pdf(file = "plots/top10_ERV_heatmap_chromosomesOnly_210620.pdf", height = 4.4, width = 3.75)
p
dev.off()

## Get counts per clade


match(c("cPa219","cPa21",
        "cPa424","cPa332",
        "HML7_repbase_HERVK11DI_MER11D","cPa492",
        "cPa568","ALVJ_Z46390",
        "SIV_NC001870","BIV_NC001413",
        "HTLV1_NC001436","BLV_NC001414",
        "cPa273","PRIMA41_repbase_PRIMA41_MER41A",
        "WDSV_NC001867","Xen1_AJ506107",
        "cPa264","SnRV_NC001724",
        "HFV_NC001736","Cer1_U15406"), ervTree$tip.label)

erv_clade.df <- tibble(erv = ervTree$tip.label,
                       clade = rep(c("Beta-like 1","Beta-like 2", "Beta", "Alpha", "Lenti", "Delta", "Gamma", "Epsilon", "SnRV-like", "Spuma"),
                                   c(120,41,17,11,10,5,56,2,17,7)))

locusID_ervAssigned_all %>% 
  left_join(erv_clade.df, by = c("erv.assigned" = "erv")) %>%
  filter(locusID %in% erv_counts.df$locusID) %>% 
  group_by(clade) %>% 
  summarise(all_erv_loci = n()) %>% 
  left_join(locusID_ervAssigned_all %>%
              filter(erv.assigned %in% top10ervs$erv.assigned) %>%
              left_join(erv_clade.df, by = c("erv.assigned" = "erv")) %>%
              group_by(clade) %>% 
              summarise(top10_erv_loci = n())) %>% 
  replace_na(list(top10_erv_loci = 0))
  

##Upset plots

upset_publish.fun <- function(df,erv,setcol,setorder){
  generic.spread <- df %>% 
    filter(erv.assigned %in% erv) %>%
    ungroup() %>% 
    select((!!as.name(setcol)), locusID, locus_present) %>% 
    spread((!!as.name(setcol)), locus_present) %>% 
    drop_na()
  
  samples_species.metadata <- loci_count.erv_assigned.allERVs.spp.summary %>%
    filter(erv.assigned %in% erv) %>%
    group_by(species) %>%
    summarise(toal_loci = sum(group_counts),
              loci_per_sample = sum(ervCounts_normalized)) %>%
    left_join(sampling_counts.species) %>%
    select(species, loci_per_sample)
  
  p <-upset(
    as.data.frame(generic.spread),
    sets = rev(setorder),#without metadata
    keep.order = TRUE,
    order.by = "freq",
    # decreasing = T,
    # nintersects = 30, #Main figure n
    nintersects = 50, #supp fig n
    mb.ratio = c(0.40,0.60),
    show.numbers = F,
    text.scale = .75,
    point.size = 1.0, 
    line.size = 0.5,
    # scale.intersections = "log2",
    # number.angles = 30,
    # set.metadata = list(data = samples_islands.metadata,
    #                     plots = list(list(type = "hist", column = "loci_per_sample", assign = 20))),
    set.metadata = list(data = samples_species.metadata,
                        plots = list(list(type = "hist", column = "loci_per_sample", assign = 20)))
    # set.metadata = list(data = samples.species_islands.metadata,
    #                     plots = list(list(type = "hist", column = "loci_per_sample", assign = 20))),
    # set.metadata = list(data = samples.hybridspecies_islands.metadata,
    #                     plots = list(list(type = "hist", column = "loci_per_sample", assign = 20))),
    # queries = list(
    # list(query = elements, params = list("hf_all", TRUE), color = "red", active = FALSE),
    # list(query = elements, params = list("fixed_all", TRUE), color = "green", active = FALSE),
    # list(query = elements, params = list("hf_1", TRUE), color = "orange", active = TRUE),
    # list(query = elements, params = list("fixed_1", TRUE), color = "blue", active = FALSE))
    )
    
    return(p)
}

## Species loci upset plot
locus_freq.df <- loci_count.erv_assigned.df %>%
  group_by(species, erv.assigned, locusID) %>% 
  summarise(locus_freq = sum(locusPresent/n())) %>% 
  mutate(locus_present = as.integer(if_else(locus_freq > 0, 1, 0)))

# erv <- locus_freq.df %>% pull(erv.assigned) %>% unique()
# for (erv in top10ervs$erv.assigned){
for (erv in n.all){ #Suppfigs
  tryCatch({
p <- upset_publish.fun(locus_freq.df,erv,"species",spplist)
# pdf(paste("plots/UpSet_", erv, "_210702.pdf", sep = ""), height = 4.6, width = 4.3)
# pdf(paste("plots/210619_", erv, "_chromosomeOnly_UpSet.pdf", sep = ""), height = 6, width = 6) #without metadata
pdf(paste("plots/suppfig_Upset/210630_", erv, "_suppFig_UpSet.pdf", sep = ""), height = 8, width = 16)
print(p)
grid.text(erv,x = 0.65, y=0.95, gp=gpar(fontsize=16))
dev.off()
  },  error=function(e){})
}

## False negative histogram
locus_freq.df %>% 
  group_by(locusID) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% #View()
  filter(c >= 24) %>% 
  pull(f) %>%  
  hist(breaks = 80, 
       xlab = "ERV frequency within species", 
       ylab = "ERV loci", 
       main = "Frequency of loci for ERVs present in all species")

locus_freq.df %>% 
  group_by(locusID) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% #View()
  filter(c >= 24) %>% 
  ggplot(aes(x = f)) +
  geom_histogram(bins = 50, color = "black", fill = "#440154FF") +
  theme_minimal() +
  theme(text = element_text(size = 6))

# locus_freq.df.summary <- locus_freq.df %>% 
#   group_by(species,locusID) %>% 
#   drop_na() %>% 
#   summarise(f = mean(locus_freq), c = sum(locus_present)) %>% 
#   filter(f != 0 && c != 0) %>% 
#   group_by(locusID) %>% 
#   summarise(ff = mean(f), cc = sum(c)) %>% 
#   filter(ff > 0.001)

locus_freq.df.summary <- locus_freq.df %>% 
  group_by(locusID) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% 
  filter(f > 0.001)

# c_f <- cut(log10(locus_freq.df.summary$f)+1, 20)
c_f <- cut(locus_freq.df.summary$f, 20)
c_c <- cut(locus_freq.df.summary$c, 24)
c_f_table <- table(c_c, c_f)

pdf(paste("plots/LocusFreqSpeciesCount_3d_210702.pdf", sep = ""), height = 8, width = 8)
# pdf(paste("plots/LocusFreqSpeciesCount_3d_chromosomeOnly_210620.pdf", sep = ""), height = 8, width = 8) 
hist3D(z=c_f_table, xlab = "Present in Species", ylab = "Frequency in Species", zlab = "log ERV loci")
image2D(z=c_f_table, border="black", xlab = "Present in Species", ylab = "Frequency in Species")
dev.off()

## Species island erv breakdown for Patric hist

p <- ggplot(data=loci_count.erv_assigned.top10.spp.island.summary %>% 
              filter(species %in% c("parvulus", "scandens", "fuliginosa","fortis")),
            aes(y = island, x = species, fill = ervCounts_normalized)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  # labs(fill = "ln(ERV loci)") +
  labs(fill = "ERV loci") +
  ggtitle("Top 10 ERV loci by species") +
  facet_wrap(~erv.assigned)

pdf("plots/210610_Top10ERV_island_species_heatmap.pdf", height = 8, width = 8) 
p
dev.off()


## Locus age histograms

afage.byspecies <- locus_freq.df %>% 
  mutate(average_allele_age = if_else(locus_freq == 1, 4,
                                      if_else(locus_freq == 0, 0,
                                              4*(-(locus_freq/(1-locus_freq))*log(locus_freq))
                                      )
  )
  )

# pdf(file = "plots/cPa75_age_210702.pdf", height = 2, width = 1.5)
pdf(file = "plots/allERVs_age_210702.pdf", height = 2, width = 1.5)
ggplot(data = afage.byspecies %>% filter(average_allele_age != 0)) +
# ggplot(data = afage.byspecies %>% filter(average_allele_age != 0, erv.assigned == "cPa75")) +
  geom_histogram(aes(x = average_allele_age), 
                 bins = 8,
                 color = "black",
                 fill = "orange") +
  labs(x = "Average allele age x*Ne", title = "Age distribution of all loci") +
  theme_minimal() +
  theme(text = element_text(size=6))
dev.off()

## Scaffold vs. Chromosome locus bias
fread("/proj/snic2020-16-24/private/ERV_Finch/Data/reference/Camarhynchus_parvulus_V1.0.fasta.fai", 
                               header = F, col.names = c("seq_name", "seq_size"), 
                               select = c(1,2)) %>%
  mutate(seq_type = if_else(grepl("chr",seq_name),"chromosome","scaffold")) %>%
  group_by(seq_type) %>% 
  summarise(cum_sum = sum(seq_size))

chrOnly.top10ervs <- c("cPa199","cPa268","cPa62","cPa75","cPa264","cPa61","cPa260","cPa443","cPa345","cPa424")

ChrScaf_loci_count <- loci_count.erv_assigned.df %>% 
  select(erv.assigned, assembly, locusID) %>% 
  group_by(erv.assigned, assembly) %>% 
  distinct() %>% 
  summarise(n = n()) %>% 
  spread(assembly, n, fill = 0) %>%
  mutate(loci_chrAsm_size = chromosome/1021297805 + 9.791463e-10,
         loci_scafAsm_size = scaffold/42718559 + 2.340903e-08)

ggplot(data = ChrScaf_loci_count, aes(x = loci_chrAsm_size, y = loci_scafAsm_size)) +
  geom_point() +
  scale_x_log10(limits = c(9.791463e-10,1e-4)) +
  scale_y_log10(limits = c(9.791463e-10,1e-4)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "log loci on chromosomes, size normalized", y = "log loci on scaffolds, size normalized")

ggplot(data = ChrScaf_loci_count,aes(x=chromosome+1, y = scaffold+1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "log loci on chromosomes + 1", y = "log loci on scaffolds + 1")

## Frequency calculations
loci_count.erv_assigned.df %>% filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
  group_by(erv.assigned) %>% 
  summarise(s = sum(locusPresent),sites = n(),freq = sum(locusPresent)/n()) %>% 
  left_join(loci_count.erv_assigned.df %>% filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
              group_by(erv.assigned) %>%
              select(erv.assigned, locusID) %>% 
              distinct() %>% 
              summarise(loci = n()))

loci_count.erv_assigned.df %>% 
  group_by(species, locusID) %>% 
  filter(erv.assigned == "cPa428") %>% 
  summarise(s = sum(locusPresent)) %>% 
  filter(s > 0) %>% 
  group_by(locusID) %>% 
  summarise(n = n()) %>% 
  group_by(n) %>% tally() #%>% 
  pull(nn) #%>% sum()

loci_present_all_spp <- loci_count.erv_assigned.df %>% 
  group_by(species, locusID) %>% 
  summarise(s = sum(locusPresent)) %>% 
  filter(s > 0) %>% 
  group_by(locusID) %>% 
  summarise(n = n()) %>%
  filter(n == 24) %>% 
  pull(locusID)

loci_present_single_spp <- loci_count.erv_assigned.df %>% 
  group_by(species, locusID) %>% 
  summarise(s = sum(locusPresent)) %>% 
  filter(s > 0) %>% 
  group_by(locusID) %>% 
  summarise(n = n()) %>%
  filter(n == 1) %>% 
  pull(locusID)

loci_count.erv_assigned.df %>% 
  filter(locusID %in% loci_present_all_spp) %>% 
  summarise(f = sum(locusPresent)/n())

loci_count.erv_assigned.df %>% 
  filter(locusID %in% loci_present_single_spp) %>% 
  group_by(species, locusID) %>%
  summarise(f = sum(locusPresent)/n()) %>% 
  filter(f != 0) %>% 
  ungroup() %>% 
  summarise(m = mean(f))

loci_count.erv_assigned.df %>% 
  filter(erv.assigned == "cPa75",
         locusPresent == 1) %>%
  mutate(incomplex = if_else(species %in% spp.list$complex, 1, 0)) %>% 
  group_by(incomplex) %>% 
  tally()

# Average frequency of top10 ERVs seggregating loci
erv_counts.df %>% 
  group_by(locusID) %>%
  summarise(locus_frequency = sum(locusPresent)/292) %>%
  left_join(locusID_ervAssigned_all) %>%
  ungroup() %>% 
  filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
  group_by(erv.assigned) %>% 
  summarise(m = mean(locus_frequency))

erv_counts.df %>% 
  select(species, island, idNumber) %>% 
  filter(species %in% spp.list$complex) %>% 
  distinct() %>% 
  summarise(n = n())

erv_counts.df %>% 
  filter(species %in% spp.list$complex) %>% 
  group_by(locusID) %>%
  summarise(locus_frequency = sum(locusPresent)/154) %>%
  left_join(locusID_ervAssigned_all) %>%
  ungroup() %>% 
  filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
  group_by(erv.assigned) %>% 
  summarise(m = mean(locus_frequency))

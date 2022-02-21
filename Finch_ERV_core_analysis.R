## Finch ERV core analysis script

library(tidyverse)
library(data.table)
library(viridis)
library(grid)
library(UpSetR)

### Data sources
#   
#   The merged dataframes are:
#   
# 1005.merged.finch.calling.csv - the merged calling matrix for retroseq, delly and CNVnator (NOTE: that delly and cnvnator are calling deletions - you still have to invert for erv calling)
# 
# 1005.merged.finch.loci.erv.assigned.csv - the list of loci with erv assignment
# 
# 
# The other files there are with more information about the retroseq loci:
#   
# 1005.all.finch.granges.1001prune.fl8clips.pad100.countfl5q1.ervoverlaps.csv - locus information about overlaps to repeatmasker calls and retrotector erv calls
# 
# 1005.all.finch.granges.1001prune.fl8clips.pad100.countfl5q1.ervsummary.csv - more information about the erv assignments to the retroseq loci
###

###
### Read in and format basic data
###

spplist <- c("noctis" , "bicolor", "olivacea","fusca","crassirostris", 
             "inornata", "parvulus",
             "pauper","heliobates",
             "pallidus","psittacula",
             "difficilis", "fuliginosa" ,"acutirostris", "BB", 
             "fortis","fortisxfuliginosa","forxfuligxscandens", "forxsca", "conirostris",
             "magnirostris","scandens", "septentrionalis","propinqua")

rs_cnv_delly.merged_call_matrix.matrix <- read_csv("1020.merged.finch.matrix.prunedERVs.csv") %>% 
  rename(locus = 1) %>% 
  rename(fortis_Daphne_15170 = BB_Daphne_15170,conirostris_Espanola_5110 = BB_Daphne_5110) %>% 
  mutate(locus = gsub("DEL[0-9]+", "DELLY", locus))

coverage.df <- read_csv("finch_292samples_STcov_reads.csv") %>% 
  mutate(samp = gsub("BB_Daphne_15170", "fortis_Daphne_15170", samp)) %>% 
  mutate(samp = gsub("BB_Daphne_5110", "conirostris_Espanola_5110", samp)) %>%
  separate(samp, into = c("species","island","idNumber"), sep = "_") %>% 
  mutate(species = factor(species, levels = spplist))


locus_ervAssigned.df <- read_csv("1020.merged.finch.loci.erv.assigned.csv") %>% 
  mutate(erv.assigned = str_extract(erv.assigned,"[a-zA-Z0-9_-]+")) %>% 
  select(locus, erv.assigned) %>% 
  mutate(locus = gsub("DEL[0-9]+", "DELLY", locus)) %>% 
  mutate(erv.assigned = gsub("HERV_","HERV-",erv.assigned))

clade_assignment.df <- read_csv("new.erv.lib.clades.csv") %>% 
  rename(erv.assigned = ERV)

sample_locus_ervAssigned.df <- rs_cnv_delly.merged_call_matrix.matrix %>% 
  gather("sampleID","count",2:ncol(.)) %>% 
  left_join(locus_ervAssigned.df) %>% 
  separate(locus, into = c("method","locus"), sep = "_") %>%
  separate(sampleID, into = c("species","island","idNumber"), sep = "_") %>%
  # filter(method == "RETROSEQ") %>% # Keep until del erv assignments are fixed
  # mutate(count_test = count) %>%
  mutate(count = if_else(method == "RETROSEQ", count,
                         if_else(count == 0, 1, 0, missing = 1))) %>%
  filter(locus != "Scaffold940:2658-10104") ## Bad Delly locus


sample_locus_ervAssigned.df %>%  filter_all(any_vars(is.na(.)))



### ERV Phylogeny

ervTree <- ape::read.nexus("camPar1_RetroAligner_AAguided_trim10gap45.pruned-2xLTRsNoSplitScaffoldERVs.FastTree2gtr.converted.Collapsed.tre")

sample_locus_ervAssigned.df %>% 
  left_join(clade_assignment.df) %>% 
  select(locus, erv.assigned, CLADE) %>% 
  distinct() %>% 
  group_by(CLADE) %>% 
  tally()

###
### General analysis
###


erv_hits.summary <- sample_locus_ervAssigned.df %>% 
  group_by(erv.assigned) %>% 
  filter(count != 0) %>% 
  tally() %>% 
  rename(erv.insertions = n) %>% 
  arrange(-erv.insertions) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label))

top10ervs <- erv_hits.summary %>% 
  filter(erv.insertions >= 10500) %>%
  arrange(-erv.insertions) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = erv.assigned))

## Summary with normalized count by sample size
sampling_counts.species <- sample_locus_ervAssigned.df %>% 
  select(species, island, idNumber) %>%
  distinct() %>% 
  group_by(species) %>% 
  summarise(birds_sampled = n())

loci_count.erv_assigned.spp.summary <- sample_locus_ervAssigned.df %>% 
  filter(count == 1) %>% 
  group_by(species, erv.assigned) %>% 
  summarise(group_counts = sum(count)) %>% 
  left_join(sampling_counts.species, by = c("species")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(erv_hits.summary$erv.assigned)))

p <- ggplot(loci_count.erv_assigned.spp.summary,
            aes(fill = (ervCounts_normalized), x = species, y = erv.assigned)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  labs(fill = "ERV insertions\nper bird") +
  ggtitle("ERV loci by species")

pdf(file = "plots/211105_All_ERV_heatmap.pdf", height = 20, width = 4.75)
p
dev.off()

### Species insertion bias
total_loci_per_sample_by_erv <- sample_locus_ervAssigned.df %>% 
  group_by(erv.assigned) %>% 
  summarise(total_insertions = sum(count), total_insertions_per_sample = total_insertions/293)

ERVinsertion_spp_bias <- loci_count.erv_assigned.spp.summary %>% 
  group_by(erv.assigned,species) %>%
  summarise(total_species_loci = sum(group_counts),
            loci_per_species_sample = sum(ervCounts_normalized)) %>%
  left_join(total_loci_per_sample_by_erv) %>% 
  mutate(normalized_insertions_per_species = loci_per_species_sample/total_insertions_per_sample)



p <-
  ggplot(ERVinsertion_spp_bias %>% filter(total_species_loci > 200)) +
  geom_boxplot(aes(x = species, y = normalized_insertions_per_species)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6))

pdf(file = "plots/211105_ERV_insertionBias_boxplot_normalized.pdf", height = 4.75, width = 4.75)
p
dev.off()

p <-
  ggplot(coverage.df) +
  geom_boxplot(aes(x = species, y = cov)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6))

pdf(file = "plots/211105_coverage_boxplot.pdf", height = 4.75, width = 4.75)
p
dev.off()

total_erv_insertions_per_species.df <- sample_locus_ervAssigned.df %>% 
  group_by(species) %>% 
  summarise(total_insertions = sum(count))

erv_proportion_per_species.df <- sample_locus_ervAssigned.df %>% 
  group_by(species, erv.assigned) %>% 
  summarise(erv_insertions = sum(count)) %>% 
  left_join(total_erv_insertions_per_species.df) %>% 
  mutate(erv_insertion_pct = erv_insertions/total_insertions)

dist_erv_proportion_per_species.mat <- erv_proportion_per_species.df %>% 
  select(species, erv.assigned, erv_insertion_pct) %>%
  spread(species, erv_insertion_pct) %>% 
  column_to_rownames(var = "erv.assigned") %>% 
  dist() 

clust_erv_proportion_per_species.mat <- hclust(dist_erv_proportion_per_species.mat)
clust_erv_proportion_per_species.mat$labels[clust_erv_proportion_per_species.mat$order]

p <- ggplot(erv_proportion_per_species.df %>%
              filter(erv.assigned %in% top10ervs$erv.assigned)) +
  geom_tile(aes(x = factor(species, levels = spplist),
                y = factor(erv.assigned, 
                           levels = (clust_erv_proportion_per_species.mat$labels[clust_erv_proportion_per_species.mat$order])),
                fill = erv_insertion_pct)) +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  ggtitle("Percentage ERV insertions of ERV group by species") +
  labs(x = "Species", y = "ERV group")

pdf(file = "plots/211105_percent_insertions_ERVgroup_species_heatmap.pdf", height = 4.75, width = 4.75)
p
dev.off()

### Top 10 ERV loci

loci_count.erv_assigned.top10.spp.summary <- sample_locus_ervAssigned.df %>% 
  filter(erv.assigned %in% top10ervs$erv.assigned) %>%
  filter(count == 1) %>% 
  group_by(species, erv.assigned) %>% 
  summarise(group_counts = sum(count)) %>% 
  left_join(sampling_counts.species, by = c("species")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(ervTree$tip.label[ervTree$tip.label %in% top10ervs$erv.assigned])))

p <- ggplot(loci_count.erv_assigned.top10.spp.summary,
            aes(fill = (ervCounts_normalized), x = species, y = erv.assigned)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  labs(fill = "ERV insertions\nper bird") +
  ggtitle("Top 10 ERV loci by species")

pdf(file = "plots/211105_top10_ERV_heatmap.pdf", height = 4.4, width = 3.75)
p
dev.off()

###
### Island analysis
###

## End goal is to look at variance of frequencies between islands where presumably
## high variance will be ERVs that have an interesting island infection history

sampling_counts.species.islands <- sample_locus_ervAssigned.df %>% 
  select(species, island, idNumber) %>%
  distinct() %>% 
  group_by(species, island) %>% 
  summarise(birds_sampled = n())

loci_count.erv_assigned.spp.island.summary <- sample_locus_ervAssigned.df %>% 
  filter(count == 1) %>% 
  group_by(species, island, erv.assigned, locus) %>% 
  summarise(group_counts = sum(count)) %>% 
  left_join(sampling_counts.species.islands, by = c("species","island")) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(ervCounts_normalized = group_counts/birds_sampled) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(erv_hits.summary$erv.assigned)))

spp_island.MinMaxVar.summary <- loci_count.erv_assigned.spp.island.summary %>% 
  group_by(species, erv.assigned, locus) %>% 
  summarise(group_count_sum = sum(group_counts),
            birds_sampled_sum = sum(birds_sampled),
            max_insertions_per_bird = max(ervCounts_normalized),
            min_insertions_per_bird = min(ervCounts_normalized),
            var_insertions_per_bird = var(ervCounts_normalized))

ggplot(spp_island.MinMaxVar.summary) +
  geom_violin(aes(x = erv.assigned, y = var_insertions_per_bird))


## Per individual sample ERV distribution
total_loci_per_individual.df <- sample_locus_ervAssigned.df %>% 
  mutate(full_sample_ID = paste(species, island, idNumber, sep = "_")) %>% 
  group_by(full_sample_ID) %>% 
  summarise(total_insertions_in_individual = sum(count))

erv_insertions_per_individual <- sample_locus_ervAssigned.df %>% 
  mutate(full_sample_ID = paste(species, island, idNumber, sep = "_")) %>% 
  group_by(full_sample_ID, erv.assigned) %>% 
  summarise(erv_insertions = sum(count)) %>% 
  left_join(total_loci_per_individual.df) %>% 
  mutate(fraction_erv_insertions_individual = erv_insertions/total_insertions_in_individual) %>% 
  separate(full_sample_ID, c("species","island","idNumber"), remove = FALSE)

p <- ggplot(erv_insertions_per_individual %>% 
              mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
              mutate(species = factor(species, levels = spplist)) %>%
              filter(erv.assigned %in% top10ervs$erv.assigned,
                     species %in% c("inornata", "noctis", "olivacea", "fortis", "scandens", "parvulus"))) +
  geom_boxplot(data = . %>% select(-species),
               aes(x = erv.assigned, y = fraction_erv_insertions_individual),
               color = "red", alpha = .5) +
  geom_boxplot(aes(x = erv.assigned, y = fraction_erv_insertions_individual),
               alpha = 0.65) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),text = element_text(size=6)) +
  facet_wrap(~species)

pdf(file = "plots/Individual_Insertion_freq_6pops_211207.pdf", height = 3.5, width = 5)
p
dev.off()

erv_insertions_per_individual_median_ratio.df <- erv_insertions_per_individual %>% 
  group_by(erv.assigned, species) %>% 
  summarise(species_insertion_median = median(fraction_erv_insertions_individual)) %>% 
  left_join(erv_insertions_per_individual %>% 
              group_by(erv.assigned) %>% 
              summarise(global_insertion_median = median(fraction_erv_insertions_individual))) %>% 
  mutate(insertion_median_ratio = species_insertion_median/global_insertion_median) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  mutate(erv.assigned = factor(erv.assigned, levels = rev(ervTree$tip.label))) 


p <- ggplot(erv_insertions_per_individual_median_ratio.df %>% 
              filter(erv.assigned %in% top10ervs$erv.assigned)) +
  geom_tile(aes(x = species, y = erv.assigned, fill = insertion_median_ratio)) +
  scale_fill_viridis(option = "D") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  labs(fill = "ERV insertion\nratio") +
  labs(title = "Median ratio of ERV insertions per individual",
       x = "ERV",
       y = "Insertion ratio")

pdf(file = "plots/Median_Insertion_Ratio_top10_colD_211116.pdf", height = 4, width = 4)
p
dev.off()

p <- ggplot(erv_insertions_per_individual_median_ratio.df %>% 
              filter(erv.assigned %in% (erv_hits.summary %>% 
                                          pull(erv.assigned) %>% 
                                          as.vector() %>% 
                                          .[1:50]))) +
  geom_tile(aes(x = species, y = erv.assigned, fill = insertion_median_ratio)) +
  scale_fill_viridis(option = "A") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  labs(fill = "ERV insertion\nratio") +
  ggtitle("ERV median insertion ratio by species")

pdf(file = "plots/211105_Median_Insertion_Ratio_top50.pdf", height = 12, width = 4)
p
dev.off()

erv_insertions_per_individual_island_median_ratio.df <- erv_insertions_per_individual %>% 
  group_by(erv.assigned, species, island) %>% 
  summarise(species_island_insertion_median = median(fraction_erv_insertions_individual)) %>% 
  left_join(erv_insertions_per_individual_median_ratio.df) %>% 
  mutate(island_species_median_ratio = species_island_insertion_median/species_insertion_median) %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = rev(ervTree$tip.label)))


p <- ggplot(erv_insertions_per_individual_island_median_ratio.df %>% 
              filter(erv.assigned %in% top10ervs$erv.assigned) %>% 
              filter(species %in% (sampling_counts.species.islands %>%
                                     group_by(species) %>%
                                     tally() %>%
                                     filter(n > 1) %>%
                                     pull(species)))) +
  # filter(species %in% c("magnirostris", "fortis", "scandens", "parvulus"))) +
  geom_tile(aes(x = island, y = erv.assigned, fill = island_species_median_ratio)) +
  geom_text(data = sampling_counts.species.islands %>% 
              filter(species %in% (sampling_counts.species.islands %>%
                                     group_by(species) %>%
                                     tally() %>%
                                     filter(n > 1) %>%
                                     pull(species))),
            # filter(species %in% c("magnirostris", "fortis", "scandens", "parvulus")),
            aes(x = island, y = "cPa183", label = paste("n=",birds_sampled,sep="")),
            color = "white",
            size = 2) +
  facet_wrap(~species, scales = 'free_x') +
  scale_fill_viridis(option = "D") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  theme(text = element_text(size = 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(fill = "ERV insertion\nratio",
       title = "ERV median insertion ratio by species between islands",
       y = "ERV") 


pdf(file = "plots/Median_Insertion_Ratio_islands_top10_211213.pdf", height = 4, width = 4)
p
dev.off()

## Species unique ERV loci

p <- sample_locus_ervAssigned.df %>%
  group_by(species, method, erv.assigned, locus) %>% 
  summarise(locus_freq = sum(count/n())) %>% 
  mutate(locus_present = as.integer(if_else(locus_freq > 0, 1, 0))) %>% 
  filter(locus_present != 0) %>% 
  group_by(method, locus) %>% 
  summarise(n_species_locus_present = n()) %>% 
  ggplot() +
  geom_histogram(aes(x = n_species_locus_present),
                 bins = 24,
                 fill = "#A3A500",
                 color = "black") +
  facet_wrap(~method, scales = "free_y")


pdf(file = "plots/211105_loci_distribution_nSpecies_method.pdf", height = 4, width = 8)
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

pdf(file = "plots/allERVs_age_211116.pdf", height = 4, width = 4)
pdf(file = "plots/top10ERVs_age_211116.pdf", height = 6, width = 6)
# ggplot(data = afage.byspecies %>% filter(average_allele_age != 0) %>% filter(erv.assigned %in% top10ervs$erv.assigned)) +
ggplot(data = afage.byspecies %>% filter(average_allele_age != 0)) +
  geom_histogram(aes(x = average_allele_age), 
                 bins = 8,
                 color = "black",
                 fill = "orange") +
  # facet_wrap(~erv.assigned, scales = "free_y") +
  labs(x = "Average allele age x*Ne", title = "Age distribution of all loci") +
  theme_minimal() +
  theme(text = element_text(size=6))
dev.off()

## False negative histogram

locus_freq.df %>% 
  group_by(locus) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% #View()
  filter(c >= 24) %>% view()
pull(f) %>%  
  hist(breaks = 80, 
       xlab = "ERV frequency within species", 
       ylab = "ERV loci", 
       main = "Frequency of loci for ERVs present in all species")

p <-locus_freq.df %>% 
  group_by(locus) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% #View()
  filter(c >= 24) %>% 
  ggplot(aes(x = f)) +
  geom_histogram(bins = 25, color = "black", fill = "#440154FF") +
  theme(text = element_text(size = 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "Insertion frequency accross all birds",
       y = "ERV loci",
       title = "Frequency of insertions of ERVs present in all species")

pdf(file = "plots/false_negative_histogram_211116.pdf", height = 4, width = 4)
p
dev.off()


locus_freq.df.summary <- locus_freq.df %>% 
  group_by(locus) %>% 
  drop_na() %>% 
  summarise(f = mean(locus_freq), c = sum(locus_present)) %>% 
  filter(f > 0.001)

## General stats
sample_locus_ervAssigned.df %>% filter(count == 1) %>% pull(locus) %>% unique() %>% length() ## Number of loci

sample_locus_ervAssigned.df %>% filter(count == 1) %>% pull(locus) %>% length() ## Number of insertions

ERV_insertion_stats.df <- tibble(erv.assigned = sample_locus_ervAssigned.df %>% 
                                   pull(erv.assigned) %>% 
                                   unique) %>% 
  left_join(clade_assignment.df) %>% 
  left_join(sample_locus_ervAssigned.df %>%  ## Loci per ERV
              filter(count==1) %>% 
              select(erv.assigned, locus) %>% 
              distinct %>%
              group_by(erv.assigned) %>% 
              tally(name = "loci")) %>% 
  left_join(sample_locus_ervAssigned.df %>%  ## Insertions per ERV
              filter(count==1) %>% 
              select(erv.assigned, count) %>% 
              group_by(erv.assigned) %>%
              summarise(insertions = sum(count))) %>% 
  left_join(sample_locus_ervAssigned.df %>%  ## species present
              filter(count==1) %>% 
              select(erv.assigned, species) %>% 
              distinct() %>% 
              group_by(erv.assigned) %>% 
              tally(name = "species_present")) %>% 
  left_join(sample_locus_ervAssigned.df %>% 
              filter(count==1) %>% 
              select(erv.assigned,locus,count) %>% 
              group_by(erv.assigned,locus) %>% 
              summarise(locus_insertion_freq = sum(count)/292) %>% 
              group_by(erv.assigned) %>% 
              summarise(locus_insertion_freq_mean = mean(locus_insertion_freq),
                        locus_insertion_freq_var = var(locus_insertion_freq)))

p <- 
  ggplot(data = ERV_insertion_stats.df) + ## Relationship between lumber of loci and number of insertions
  geom_point(aes(x = log(loci), y = log(insertions))) +
  labs(title = "ERV loci vs. identifications",
       x = "log ERV loci",
       y = "log ERV identifications")
pdf(file = "loci_vs_insertions.pdf", width = 4, height = 4)
p
dev.off()

ERV_insertion_stats.df %>% 
  filter(grepl("cPa", erv.assigned)) %>% 
  summarise(max_insertions = max(insertions),
            min_insertions = min(insertions),
            mean_insertions = mean(insertions),
            median_insertions = median(insertions),
            std_deviation_insertions = sd(insertions),
            max_loci = max(loci),
            min_loci = min(loci),
            mean_loci = mean(loci),
            median_loci = median(loci),
            std_deviation_loci = sd(loci),
            mean_freq = mean(locus_insertion_freq_mean),
            sd_freq = sd(locus_insertion_freq_mean)) %>% View()


ERV_insertion_stats.df %>%
  select(erv.assigned, species_present) %>% 
  filter(grepl("cPa",erv.assigned)) %>% 
  group_by(species_present) %>% 
  tally()

ggplot(data = locus_freq.df %>% 
         filter(locus_freq > 0, erv.assigned %in% top10ervs$erv.assigned) %>% 
         left_join(sampling_counts.species) %>% 
         mutate(species = paste(species, birds_sampled, sep = " "))) +
  geom_violin(aes(x = erv.assigned, y = locus_freq, fill = species)) +
  facet_wrap(~erv.assigned, scales = "free_x")

locus_freq.df %>% ## Number of loci fixed in all species
  filter(locus_freq >= 1, grepl("cPa", erv.assigned)) %>% 
  group_by(erv.assigned,locus) %>%
  tally() %>% 
  arrange(-n) %>% 
  filter(n == 24) %>% 
  group_by(erv.assigned) %>% 
  tally() %>% 
  arrange(-n) %>% view()

locus_freq.df %>% ## Number of loci fixed in at least one species
  filter(locus_freq >= 1, grepl("cPa", erv.assigned)) %>% 
  select(erv.assigned, locus) %>% 
  group_by(erv.assigned) %>% 
  distinct() %>% 
  tally() %>%
  view()

fixed_loci_ERV_species.df <- locus_freq.df %>% ## Number of loci fixed by ERV and species
  filter(locus_freq >= 1, grepl("cPa", erv.assigned)) %>% 
  group_by(erv.assigned,locus, species) %>%
  tally() %>% 
  group_by(erv.assigned, species) %>% 
  summarise(loci_fixed = sum(n))

coverage_freq_insetion_ratio.comapre.df <- sample_locus_ervAssigned.df %>% 
  filter(grepl("cPa", erv.assigned)) %>% 
  group_by(species, island, idNumber) %>% 
  summarise(cPa_freq = sum(count)/n()) %>% 
  left_join(erv_insertions_per_individual %>% ## Normalized coverage vs. individual insertion ratio
              filter(grepl("cPa", erv.assigned)) %>% 
              group_by(full_sample_ID, species, island, idNumber) %>% 
              summarise(median_fraction_insertion_individual = median(fraction_erv_insertions_individual))) %>% 
  left_join(coverage.df) 

coverage_freq_insetion_ratio.comapre.df %>% ggplot() + 
  geom_point(aes(x = cov, y = cPa_freq)) +
  labs(x= "Individual sequencing coverage", y = "Median ERV frequency")

lf <- lm(cov ~ cPa_freq, coverage_freq_insetion_ratio.comapre.df)
lf <- lm(cov ~ median_fraction_insertion_individual, coverage_freq_insetion_ratio.comapre.df)

## Fig 4 stats
fig4.aov <- aov(fraction_erv_insertions_individual ~ species, erv_insertions_per_individual %>% 
                  mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
                  mutate(species = factor(species, levels = spplist)) %>%
                  filter(erv.assigned == "cPa260",
                         species %in% c("inornata", "noctis", "olivacea", "fortis", "scandens", "parvulus")))
summary(fig4.aov)
TukeyHSD(fig4.aov)

cPa363.aov <- aov(fraction_erv_insertions_individual ~ species, erv_insertions_per_individual %>% 
                    mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
                    mutate(species = factor(species, levels = spplist)) %>%
                    filter(erv.assigned == "cPa363"))
summary(cPa363.aov)
cPa363.aov.tukeyHSD <- TukeyHSD(cPa363.aov)
cPa363.aov.tukeyHSD$species %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "comparison") %>%
  rename(p_adj = `p adj`) %>% 
  arrange(p_adj) %>% 
  write_tsv(file = "cPa363_tukeyHSD.tsv")

t.test(x = erv_insertions_per_individual %>% 
         mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
         mutate(species = factor(species, levels = spplist)) %>%
         filter(erv.assigned == "cPa363",
                species == "olivacea") %>% pull(fraction_erv_insertions_individual),
       y =  erv_insertions_per_individual %>% 
         mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
         mutate(species = factor(species, levels = spplist)) %>%
         filter(erv.assigned == "cPa363",
                species != "olivacea") %>% pull(fraction_erv_insertions_individual))

fig4.aov <- aov(fraction_erv_insertions_individual ~ species, erv_insertions_per_individual %>% 
                  mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
                  mutate(species = factor(species, levels = spplist)) %>%
                  filter(erv.assigned == "cPa363",
                         species != "olivacea"))

erv_insertions_per_individual %>% 
  mutate(erv.assigned = factor(erv.assigned, levels = ervTree$tip.label)) %>% 
  mutate(species = factor(species, levels = spplist)) %>%
  filter(erv.assigned == "cPa363",
         species != "olivacea") %>% pull(fraction_erv_insertions_individual) %>% mean()

locus_freq.df %>%
  filter(species == "olivacea", erv.assigned == "cPa363") %>% 
  filter(locus_freq > 0) %>% 
  summarise(m = mean(locus_freq), sd = sd(locus_freq))

## ERV density analysis
genome_assembly_index.df <- read_table("Camarhynchus_parvulus_V1.0.fasta.fai", 
                                       col_names = c("chromosome", "size", "x", "y", "z")) %>% 
  select(chromosome, size) %>% 
  mutate(chromosome = gsub("Scaffold[0-9]+", "Scaffold", chromosome)) %>% 
  group_by(chromosome) %>% 
  summarise(size = sum(size))

total_loci_per_individual_per_erv.df <- sample_locus_ervAssigned.df %>% 
  group_by(species, 
           island,
           idNumber,
           erv.assigned) %>% 
  summarise(sample_total_genome = sum(count))

erv_density_per_individual_per_chromosome.df <- sample_locus_ervAssigned.df %>% 
  separate(locus, into = c("chromosome", "coordinate"), sep = ":") %>% 
  select(chromosome, species, island, idNumber, count, erv.assigned) %>% 
  mutate(chromosome = gsub("Scaffold[0-9]+", "Scaffold", chromosome)) %>% 
  group_by(chromosome, species, island, idNumber, erv.assigned) %>% 
  summarise(sample_total_chromosome = sum(count)) %>% 
  left_join(total_loci_per_individual_per_erv.df) %>% 
  left_join(genome_assembly_index.df)

ggplot(erv_density_per_individual_per_chromosome.df %>% 
         filter(sample_total_genome > 10, sample_total_chromosome >= 1) %>% 
         filter(erv.assigned %in% top10ervs$erv.assigned)) +
  geom_point(aes(x = sample_total_genome, 
                 y = sample_total_chromosome,
                 color = erv.assigned)) 

genome_assembly_index.df %>% summarise(genome_size = sum(size)) # 1064016364

erv_density_per_individual_per_chromosome.df <- erv_density_per_individual_per_chromosome.df %>% 
  mutate(ERV_chromosome_density_per_MB = (sample_total_chromosome/size)*1000000) %>% 
  mutate(ERV_genome_density_per_MB = (sample_total_genome/1064016364) * 1000000) %>% 
  mutate(ERV_density_deviation_from_genome_average = ERV_chromosome_density_per_MB/ERV_genome_density_per_MB)


ggplot(erv_density_per_individual_per_chromosome.df %>% 
         filter(erv.assigned %in% top10ervs$erv.assigned))+
  geom_point(aes(x = ERV_genome_density_per_MB, 
                 y = ERV_chromosome_density_per_MB, 
                 color = erv.assigned)) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~chromosome, scales = "free_y")

## Investigation of distribution of specific ERVs shown to deviate above

density_map.df <- sample_locus_ervAssigned.df %>% 
  separate(locus, into = c("chromosome", "coordinate"), sep = ":") %>% 
  separate(coordinate, into = c("start", "end"), sep = "-") %>% 
  filter(erv.assigned %in% top10ervs$erv.assigned,
         grepl("chr", chromosome)) %>% 
  select(chromosome, start, end, erv.assigned) %>% 
  distinct()

density_map.df <- density_map.df %>% 
  mutate(chromosome = factor(chromosome, levels = (genome_assembly_index.df %>% 
                                                     arrange(-size) %>% 
                                                     pull(chromosome))),
         start = as.numeric(start),
         end = as.numeric(end)) %>% 
  left_join(genome_assembly_index.df)

ggplot(density_map.df) +
  geom_vline(aes(x = chromosome, y = size), alpha = 0.25) +
  geom_point(aes(x = chromosome, y = start, color = erv.assigned)) +
  facet_wrap(~erv.assigned)

ggplot(density_map.df) +
  geom_density(aes(x = start, fill = erv.assigned), alpha = 0.25) +
  facet_wrap(~chromosome, scales = "free")

## Density of all ERVs and cPa452 over chromosomes between 3 groups of interest

locus_freq.df.summary %>% 
  group_by() %>% 
  separate(locus, into = c("chromosome", "coordinate"), sep = ":") %>% 
  separate(coordinate, into = c("start", "end"), sep = "-") %>% 
  mutate(chromosome = factor(chromosome, levels = (genome_assembly_index.df %>% 
                                                     arrange(-size) %>% 
                                                     pull(chromosome))),
         start = as.numeric(start),
         end = as.numeric(end)) %>% 
  drop_na() %>% 
  distinct() %>% 
  ggplot() +
  geom_histogram(aes(x = start), binwidth = 1000000, center = 500000, fill = "gray70", color = "black") +
  facet_wrap(~chromosome, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  labs(x = "Genomic position", y = "ERV loci per Mb")

locus_freq.df.summary %>% 
  group_by() %>% 
  separate(locus, into = c("chromosome", "coordinate"), sep = ":") %>% 
  separate(coordinate, into = c("start", "end"), sep = "-") %>% 
  drop_na() %>% 
  distinct() %>% 
  group_by(chromosome) %>% 
  tally() %>% 
  left_join(genome_assembly_index.df) %>% 
  mutate(loci_per_Mb = (n/size) * 1000000) %>% 
  mutate(chromosome = factor(chromosome, levels = (genome_assembly_index.df %>% 
                                                     arrange(-size) %>% 
                                                     pull(chromosome)))) %>% 
  drop_na() %>% 
  write_tsv(file = "SupplementaryTable1_ERV_DensityPerChrom_220221")

erv_density_per_individual_per_chromosome.df %>% 
  mutate(spp_group = case_when(
    species %in% c("heliobates","pallidus","parvulus","pauper","psittacula") ~ "Tree",
    species %in% c("olivacea","fusca") ~ "Warbler",
    species %in% c("acutirostris","conirostris","difficilis","fortis","fuliginosa","magnirostris",
                   "propinqua","scandens","septentrionalis") ~ "Ground",
    TRUE ~ "Other")
  ) %>% 
  filter(erv.assigned == "cPa452", spp_group != "Other") %>% 
  mutate(chromosome = factor(chromosome, levels = (genome_assembly_index.df %>% 
                                                     arrange(-size) %>% 
                                                     pull(chromosome)))) %>% 
  ggplot()+
  geom_point(aes(x = ERV_genome_density_per_MB, 
                 y = ERV_chromosome_density_per_MB, 
                 color = spp_group), alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~chromosome, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),text = element_text(size=6)) +
  labs(x = "Genome wide ERV insertions per Mb", y = "ERV insertions per Mb by chromosome")
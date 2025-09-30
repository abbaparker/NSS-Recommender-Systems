#set working directory where NSS excel file is saved
setwd("")
library(tidyverse)
library(readxl)
library(corrplot)

#Read in data
d1 <- read_excel("df_species_with_100yrtime_increments_Sept2025Abba.xlsx")
colnames(d1)

# counts per region
d1 %>% dplyr::count(Region)

#For this run, remove Poland and Estonia (different datasets)
d2 <- d1 %>% filter(!(Region=="Central & Eastern Europe"))

#Species with _
d1$GBIF_species <- gsub(" ", "_", d1$GBIF_species)

# Create List of Site Metadata
assemblage <- d1 %>% dplyr::select("DB_Assemblage_ID", "Site_name", "Settlement_modern_name",
                            "Country", "Region", "time_bins2",
                            "Rural_urban_or_neither",
                            "Historic_Urban_Centre", "DB_sieved",
                            "start_date_CE_final_num", "end_date_CE_final_num", "Earliest_relevant_urban_date_in_NSS",    
                            "Decimal_Latitude", "Decimal_Longitude","Time.mid") 

fassemblage <-  unique(assemblage)  
write.csv(fassemblage, "Site_Metadata_FamilyLevel_noCEE_190925.csv")
assemblage <-  unique(assemblage)
write.csv(assemblage, "Site_Metadata_SpeciesLevel_noCEE_190925.csv")
# this will vary for species data (different sites excluded as top/bottom 2.5% of counts)
    # so re-run

sitec <- assemblage %>% count(DB_Assemblage_ID) #Check-- is each assemblage ID unique?
assemblage %>% count(DB_sieved) #how many sieved assemblages
assemblage %>% count(time_bins2) # coutn by time bin
assemblage %>% count(Historic_Urban_Centre) # count urban/nonurban assembalges
checkf <- assemblage %>% filter(Historic_Urban_Centre=="NA") #what is missing this coding?
#write.csv(checkf,"NAurbancenters.csv")

# Create spreadsheet of taxon metadata

# For family data
# find families with over 3 occurrences
dfam <- d2 %>% dplyr::select(GBIF_family, GBIF_order)
occ_count <- dfam %>% dplyr::count(GBIF_family)
fam3plus <- occ_count %>% filter(n>2) %>% filter(!is.na(GBIF_family))
ftaxa <- fam3plus$GBIF_family
ftaxa <- ftaxa[!(ftaxa=="NA")]  # remove NA's
ftaxa #results in 52 families to use

# Save family trait data
data_families <- d2 %>% select("GBIF_family",
                               "GBIF_order","Trait_Habitat",   "Trait_LifeHistory", 
                               "Trait_Temp_Max", "Trait_Temp_Min", 
                               "Known_Aquaculture_R_Hoffmann", 
                               "Known_Major_Trade_J_Barrett")
data_families <- unique(data_families)
write.csv(data_families, "TraitData_Families_NSS_19Feb25.csv")

#for species data
dsp <- d1 %>% dplyr::select(GBIF_genus, GBIF_family, GBIF_order, GBIF_species)
occ_count <- dsp %>% dplyr::count(GBIF_species)  
taxa3plus <- occ_count %>% filter(n>2)  # species with 3 or more occurrences
taxa <- taxa3plus$GBIF_species
taxa #113 species, 116 with CEE

#save species trait data
data_taxa_all <- d1 %>% dplyr::select("GBIF_species", "GBIF_genus", "GBIF_family",
                               "GBIF_order",  "GBIF_class", "GBIF_phylum" ,"GBIF_kingdom", "Trait_Habitat",         
                                 "Trait_LifeHistory",
                               "Trait_Temp_Max", "Trait_Temp_Min",
                               "Known_Aquaculture_R_Hoffmann", "Known_Major_Trade_J_Barrett")
data_species <- unique(data_taxa_all)
write.csv(data_species, "TraitData_Species_NSS_19Sep25.csv")

# Build site/occurrence matrices ####

#choose one of family or species level, read above saved files 
straits <- read.csv("TraitData_Species_NSS_19Sep25.csv")
taxa <- taxa  #species
#taxa <- ftaxa #families

# select site names
sites <- read.csv("Site_Metadata_SpeciesLevel_noCEE_190925.csv", row.names = 1) #or SpeciesLevel
site <- unique(sites$DB_Assemblage_ID) 

#make NISP numeric
d1$NISP <- as.numeric(d1$NISP)

#Build matrices
Occ_matrix <- as.data.frame(matrix(nrow=length(site), ncol=length(taxa))) # presence/absence by site
Nisp_matrix <- as.data.frame(matrix(nrow=length(site), ncol=length(taxa))) # NISP of each taxon for each site
colnames(Occ_matrix) <- taxa
colnames(Nisp_matrix) <- taxa
rownames(Occ_matrix) <- site
rownames(Nisp_matrix) <- site

# Total matrix for all unique taxon names

# For species run
for (q in 1:nrow(sites)){  #row for each locality
  pt <- d1 %>% filter(DB_Assemblage_ID==sites[q,1]) 
  for (r in 1:length(taxa)){  #column for each taxon
    # find instances of that taxon from site ID
    records <- pt %>% filter(GBIF_species==taxa[r])
    if (nrow(records) >0) {# where there are records 
      Occ_matrix[q,r] <- 1
      Nisp_matrix[q,r] <- sum(records$NISP)
    }
    else {
      Occ_matrix[q,r] <- 0
      Nisp_matrix[q,r] <- 0
    }
  }
}
write.csv(Occ_matrix, "Species_PA_noCEE_19Sep25.csv")
write.csv(Nisp_matrix, "Species_NISP_noCEE_19Sep25.csv")

#For family run
for (q in 1:nrow(sites)){  #row for each locality
  pt <- d1 %>% filter(DB_Assemblage_ID==sites[q,1]) 
  for (r in 1:length(taxa)){  #column for each taxon
    # find instances of that taxon from site ID
    records <- pt %>% filter(GBIF_family==taxa[r])# %>% filter(ID %in% pt$ID)
    if (nrow(records) >0) {# where there are records 
      Occ_matrix[q,r] <- 1
      Nisp_matrix[q,r] <- sum(records$NISP)
    }
    else {
      Occ_matrix[q,r] <- 0
      Nisp_matrix[q,r] <- 0
    }
  }
}

write.csv(Occ_matrix, "Family_PA_CEE_18Feb25.csv")
write.csv(Nisp_matrix, "Family_NISP_CEE_18Feb25.csv")


#Histogram plot-- how many taxa per site?
data_all <- read.csv("Species_PA_noCEE_19Sep25.csv", header = TRUE)
Ntax <- rowSums(data_all[2:ncol(data_all)], na.rm = T)
ggplot() + geom_histogram(aes(x=Ntax))


# For evaluation, find NISP for sieved sites only ####
colnames(d1)
sieveS <- unique(d1$DB_Assemblage_ID[which(d1$DB_sieved==T)]) #list of sites that have been sieved
length(sieveS)/length(site)  #fraction of sites that are sieved
data_nisp <- read.csv("Species_NISP_noCEE_19Sep25.csv", header = TRUE)
colnames(data_nisp)
sieve_nisp <- as.data.frame(matrix(nrow=length(site), ncol=length(taxa)))
colnames(sieve_nisp) <- taxa
rownames(sieve_nisp) <- site

for (g in 1:length(sieveS)){  # for all sieved sites
  rnum <- which(data_nisp$X == sieveS[g]) #which row of total matrix is the site
  Ndat <- data_nisp[data_nisp$X == sieveS[g],2:100]  #data from the site
  Nsum <- sum(Ndat)  #total NISP for the sieved site
  frac <- Ndat/Nsum #fraction of NISP for each taxon
  sieve_nisp[rnum,] <- frac  #save these fractions into the row
}
write.csv(sieve_nisp, "Sieved_Species_NISPfraction_noCEE_19Sep25.csv")

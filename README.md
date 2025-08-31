Code accompanying the manuscript "2000-year fish bone record reveals transition to commercial fisheries 
during climatic change over the past two millennia" (Buss et al. 2025).

Contents of R scripts:
1) Input NSS-DB data, save as assemblage (row)/taxon (column) matrices indicating databased occurrences 
2) Test for similarity between sites (labelled as historic urban centres or not) in fish bone occurrences, relative to null
3) Test for similarity between species (labelled by life histories) in their occurrence patterns, relative to null
4) Run Recommender Systems using matrix factorization, including parameter testing.

The outputs of script 4 are evaluation of performance with different parameter combinations, as well as "REC" files 
with the estimated propensity of each taxon to occur at each sampled assemblage, "Xmat" files showing the strength of
the association of each sampled assemblage to the identified fishery types (termed V1-V4 in the manuscript), and "Ymat"
files showing the association of each taxon with those same fishery types. 

Any questions? Please get in touch with abigail.parker@helsinki.fi

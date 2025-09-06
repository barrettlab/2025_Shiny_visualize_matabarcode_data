# 2025_Shiny_visualize_matabarcode_data

The code includes an R shiny app to plot fungal ITS2 metabarcoding data from a phyloseq object (.rds)

The code allows users to:

- View taxa at different rank levels (Phylum, Class, Order, etc.) interactively

- Hover over plot to see pop-up of full taxonomic information for each OTU

- Control number of columns in faceted plots based on sample_data() of phyloseq object

- Save plots as pdf

Two example datasets (.rds file extension)

#### In R, how to save and read in .rds file of phyloseq object
```{r{
# Save the object, where ps_rel is the phyloseq object in R
saveRDS(ps_rel, file = "Corallorhiza_striata_ps_rel.rds")

# Load the object
ps_rel <- readRDS("Corallorhiza_striata_ps_rel.rds")
```

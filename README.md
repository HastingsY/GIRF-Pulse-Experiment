# GIRF Pulse Experiment
R code for Hastings, Y. D. (2022). Green Infrastructure Microbial Community Response to Simulated Pulse Precipitation Events in the Semi-Arid Western United States (Master's thesis, The University of Utah). This study was supported by a grant from the US National Science Foundation (DEB 2006308).

R code for and Hastings, Y. D., et al. Green Infrastructure Microbial Community Response to Simulated Pulse Precipitation Events in the Semi-Arid Western United States. In preparation.

Abstract: Nutrient retention in urban stormwater green infrastructure (SGI) of water-limited biomes is not well quantified, especially when stormwater inputs are scarce. We examined the role of plant diversity and physiochemistry as drivers of microbial community physiology and soil N pools and fluxes in bioswales subjected to simulated precipitation and a montane meadow experiencing natural rainfall within a semi-arid region during drought. Precipitation generally elevated soil moisture and pH, stimulated ecoenzyme activity, and increased the concentration of organic matter, proteins, and N pools in both bioswale and meadow soils; but the magnitude of change differed between events. Microbial community growth was static and N assimilation into biomass was limited across precipitation events. Unvegetated SGI plots had greater soil moisture, yet effects of plant diversity treatments on microbial C:N ratios, organic matter content, and N pools were inconsistent. Differences in soil N concentrations in bioswales and the meadow were most directly correlated to changes in organic matter content mediated by ecoenzyme expression and the balance of C, N, and P resources available to microbial communities. Our results add to growing evidence that ecological function of SGI is comparable to neighboring natural vegetated systems, particularly when soil media and water availability are similar.

The file and R code structure is as follows: 

1. Data - Contains all data used for the analysis
2. Results - Contains all figures, RMANOVA, and Piecewise Structural Equation Modeling results. 
3. renv - R environment used for project
4. EEA_Vector_Analysis.R - R code used to analyze coenzyme (EEA) responses, including RMANOVA to look for significant differences in EEA response to simulated pulse events and Vector Analysis to determine the nutrient resource acquisition.
5. Gravimetric_soil_moisture_pH.R - R code used for RMANOVA of gravimetric soil moisture and pH responses to simulated pulse events.
6. MicrobialBiomass_EEA.Rproj - Downloaded R project
7. Microbial_biomass.R - R code used for RMANOVA of microbial biomass carbon, nitrogen, and C:N responses to simulated pulse events.
8. OM_protien_N_pools_fluxes.R - R code used for RMANOVA of organic matter content, proteins, and N pools and fluxes responses to simulated pulse events.
9. PSEM_final.R - R code used for Pearson Correlation and Piecewise Structural Equation Modeling.
10. Rclimate.R - R code used to obtain summary statistics of climate data from GIRF and TM climate and soil sensors.



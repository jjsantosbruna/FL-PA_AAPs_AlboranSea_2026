Scripts for Aerobic anoxygenic phototrophic bacteria correlate with picophytoplankton across the Atlantic Ocean but show unique vertical bioenergetics
A manuscript published in Limnology & Oceanography doi: https://doi.org/10.1002/lno.12682

Scripts by Carlota R. Gazulla

Data and scripts are organized as follows:

Data
A main file Metadata Gazulla et al 2024 L&O.xlsx summarizes all the metadata. The .csv files that should be used in the scripts are listed below.

SNPP_VIIRS.20190301_20190331.L3m.MO.CHL.chlor_a.9km.nc is the near-surface chlorophyll data from March 2019, which is represented in Figure 1. It was downloaded from https://oceancolor.gsfc.nasa.gov/l3/.
expedition_data.csv contains basic data from the Poseidon Expedition (stations, depths, coordinates...), cell abundances from flow cytometry and pigment concentrations from the HPLC.
light_data.csv contains each station's photosynthetic active radiation (PAR) and diffuse attenuation coefficients at 490 nm (Kd). The data has been downloaded from https://oceancolor.gsfc.nasa.gov/l3/.
integrated_data.csv contains the integrated cellular daily energy yield at each station.
Scripts
0_Bioenergetics.R contain light and bioenergetics calculations. It should be executed before the following scripts.
Figure1.R is to generate Figure 1.
Figure2.R to generate panels B (right) and C (right) in Figure 2.
Figure3.R to generate Figure 3.
Figure4.R to generate Figure 4.
Supplementary_Figures.R to generate Supplementary Figures S7, S8 and S9.

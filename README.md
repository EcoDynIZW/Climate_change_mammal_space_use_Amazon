# Rocha and Sollmann (2023) *AnimCOnserv*

> Rocha DG and **Sollmann R** (2023): Habitat use patterns suggest that climate-driven vegetation changes will negatively impact mammal communities in the Amazon. *Anim Conserv*. DOI: [10.1111/acv.12853](https://doi.org/10.1111/acv.12853)

This repository contains the R scripts and data to implement the two main community
occupancy models for mammals in the southern Brazilian Amazon described in Rocha and
Sollmann (2023). Habitat use patterns suggest that climate-driven vegetation changes will negatively impact mammal communities in the Amazon. Animal Conservation, 
https://doi.org/10.1111/acv.12853
The material is also available on Dryad under Sollmann and Rocha (2023): Data and model code for: Habitat use patterns suggest that climate-driven vegetation changes will negatively impact mammal communities in the Amazon (ACV), Dryad, Dataset, https://doi.org/10.25338/B84060 
Specifically, the data are camera-trap data from 5 surveys in 4 protected areas in the southern Brazilian Amazon, summarized as species detection/non-detection data, along with camera-specific covariates and effort information.


## Description of the data and file structure

The file 'Rocha Sollmann data.R' contains a list with all data necessary to fit the two community occupancy models described in the main text of Rocha and Sollmann 2023, Habitat use patterns suggest that climate-driven vegetation changes will negatively impact mammal communities in the Amazon, Animal Conservation, https://doi.org/10.1111/acv.12853.

The list has the following elements:
$ Y: array, location by occasion by species-survey, with entries of 1 denoting detection of a species, 0 denoting non-detection, and NA for location-occasions combinations that were not sampled. Data from the 5 surveys are stacked.
$ covmat: array, survey by location by covariate, with location-specific covariate values (P.Savanna = proportion savanna habitat in 100-m buffer; Hab.type = 0/closed or 1/open; dRio = scaled distance to nearest river; d.PA = scaled distance to protected area border; pWet = proportion of survey in the wet season; sav.sc.ht = proportion of savanna habitat in 100-m buffer, scaled by immediate habitat type; mata = immediate habitat is continuous forest; ripa = immediate habitat is riparian forest; cerr = immediate habitat is savanna; floresta = immediate habitat is forest; sav.sc.ht.01 = not relevant; savG.sc = scaled proportion of savanna habitat in 100-m buffer)
$ pWET.occ: array, survey by location by occasion, with proportion of sampling occasion in the wet season
$ species: vector, numeric code for species identity
$ survey: vector, numeric code for survey identity
$ region: vector, numeric code for region (north = 1 or south = 2)
$ K: vector: number of sampling occasions in each survey
$ Jmax: scalar, maximum number of locations in any survey
$ nsite: vector, number of locations in each survey
$ n.survey: scalar, number of surveys
$ effort.sc: array, survey by location by occasion, with scaled sampling effort
$ first: array, survey by location, with first occasion with sampling effort
$ last: array, survey by location, with last occasion with sampling effort
$ n.spec: scalar, number of species
$ n.sp.surv: scalar, number of species-survey combinations (not equal to number of surveys times number of species because some species were not detected in all surveys)
$ sp.names.cor: vector with Brazilian species names, same order as in data set (English names are provided in the R code to run models)


In addition, the following code files are available:

R_code_to_run_model_1.R: R code to load data, prep it for analysis and implement the first model described in the publication, estimating the effect of proportion of savanna habitat (linear and squared) separately by region on species occupancy probability; note: requires Nimble_custom_function.R, and Nimble_model_savQ_27Oct2021.R to be located in the working directory

R_code_to_run_model_2.R: R code to load data, prep it for analysis and implement the second model described in the publication, estimating the effect of proportion of savanna habitat by immediate habitat type and species group on species occupancy probability; note: requires Nimble_custom_function.R, and Nimble_model_sav-geral-intHT_group.R to be located in the working directory

Nimble_model_savQ_27Oct2021.R: Nimble model code for model 1 (see above)

Nimble_model_sav-geral-intHT_group.R: Nimble model code for model 2 (see above)

Nimble_custom_function.R: Some custom distribution functions used in Nimble model code (vectorizing over occasions, integrating out occupancy state)


## Sharing/Access information

NA


## Code/Software

Data can be loaded into R and analyzed using the R package Nimble. Please see publication for version information. 
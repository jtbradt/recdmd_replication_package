# Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models

### Author: Jacob Bradt

The code and data contained in this directory replicate 'Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models.'

> **Abstract:** Conventional applications of recreation demand models likely suffer from two standard challenges with demand estimation, namely omitted variables bias and measurement error. Idiosyncratic prices in the form of individual-level travel costs can exacerbate these two challenges: the potential for non-random selection into travel costs through residential sorting and the difficulty of observing individual-level travel costs both work to bias traditional model estimates. I demonstrate the magnitude of this potential bias in conventional estimates of recreation demand models. I provide a relatively simple instrumental variables approach to address these two empirical challenges that substantially outperforms traditional estimates in numerical simulations. Replicating English et al. (2018), I find that accounting for potential selection into travel costs and measurement error through the instrumental variables approach decreases estimates of the welfare costs of the 2010 Deepwater Horizon oil spill by 12 percent.

The code in this replication package constructs the necessary simulated and empirical datasets used in the simulation and English et al (2018) empirical exercises contained in the paper.  Subsequent code uses these datasets to estimate the main results of the paper, generating 11 of 12 figures and 7 of 7 tables contained within the main text and online appendices.  

All code is written in R and Julia. I have tested this replication package using both Windows and Mac operating systems.

Computation of the full set of simulations reported in the text is time-intensive.  As a result, these replication files reduce the number of Monte Carlo simulations from 1000 to 100 for each set of simulated results.  Similarly, standard errors reported in the English et al. (2018) empirical exercise in the main text require a bootstrap procedure, which is computationally intensive.  I also reduce the number of bootstrap iterates from 200 to 20 in these empirical results.  Users interested in obtaining the full results from the paper are welcome to adjust the number of simulations and bootstrap iterations. Note that the versions of all figures and tables contained within this replication package (stored in the `/output/` subfolder) correspond to the final versions contained within the paper and online appendices.  Executing this replication package with the reduced number of simulations and bootstrap iterations will overwrite these figures and tables. 

With these reductions in computational burden, the replication package should take around 10 hours when run on a machine with 192GB of memory and 16 cores.

## Data Availability and Provenance Statements

-   [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

-   [x] All data **can be made** publicly available.

## Data Description

The data in this repository are organized into two separate folders, one each for the simulation results and the English et al. (2018) empirical replication exercise, respectively.  The intermediate and final simulated data generated as part of the Monte Carlo analysis are stored in [`/data/simulations/`](/data/simulations/).  The raw, intermediate, and final data used as part of the English et al. (2018) empirical exercise are stored in [`/data/dwh_replication/`](/data/dwh_replication/).

### [Simulation Data](/data/simulations/)

The data stored in the [Simulation Data](/data/simulations/) subfolder are generated directly by Julia code outlined in a separate section below.  As a result, these files store intermediate datasets and final estimates from the simulation exercises.  Additional details on specific components of this subfolder are discussed below.

1. [sim_data_sim1to6.csv](/data/simulations/sim_data_sim1to6.csv): contains single simulated data set for each of the data generating processes defined by Simulations 1-6 in the text.

2. [est_wtp_sim1.csv](/data/simulations/est_wtp_sim1.csv): contains estimated willingness-to-pay (WTP) statistics across all simulated datasets for the data generating process defined by Simulation 1 in the text.  Includes WTP estimates based on each of the three alternative estimators discussed in the text---the baseline multinomial logit estimator, the baseline multinomial logit estimator with origin group fixed effects, and two-stage control function estimator.  Similarly defined CSV files contain analogous WTP estimates for Simulations 2 through 6.

3. [est_wtp_rc1.csv](/data/simulations/est_wtp_rc1.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different degrees of within-recreator origin group correlation. The results from this robustness check are reported in Appendix Figure A1.

4. [est_wtp_rc2.csv](/data/simulations/est_wtp_rc2.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different recreator origin group sizes. The results from this robustness check are reported in Appendix Figure A2.

5. [est_wtp_rc3.csv](/data/simulations/est_wtp_rc3.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different forms of nonlinearity in the unobservable, idiosyncratic preference term. The results from this robustness check are reported in Appendix Figure A3.

6. [est_wtp_rc4.csv](/data/simulations/est_wtp_rc4.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different forms of nonlinearity in the cost instrument. The results from this robustness check are reported in Appendix Figure A4.

7. [est_wtp_rc_weak.csv](/data/simulations/est_wtp_rc_weak.csv) and [est_wtp_rc_weak_fstat.csv](/data/simulations/est_wtp_rc_weak_fstat.csv): contain results from the robustness check that implements the Simulation 2 data generating process with different degrees of instrument relevance. The results from this robustness check are reported in Appendix Figure A5.

### [Empirical Exercise Data](/data/dwh_replication/)

The data stored in the this subfolder are categorized into three separate types: (1) data taken directly from English et al. (2018), (2) supplemental data used to generate travel cost intstruments, and (3) constructed files containing intermediate data used to generate parameter estimates and standard errors.

#### English et al. (2018) source data
I obtain the source data on the Deepwater Horizon oil spill and Gulf Coast beach visitation necessary for the English et al. (2018) empirical exercise from the National Oceanographic and Atmospheric Administration's [Natural Resource Damage Assessment (NRDA) public repository](https://www.diver.orr.noaa.gov/deepwater-horizon-nrda-data).  These data are made available to the public and as stated above I have legitimate access and permission to use and report these data as part of this replication package.  I only include those source datasets from English et al. (2018) which I use directly in the production of the results reported in this paper.  These datasets include the following (see the main text or English et al. (2018) for additional information):

1. [origins_local.csv](/data/dwh_replication/origins_local.csv) and [origins_national.csv](/data/dwh_replication/origins_national.csv): contain respondent origin location information for all respondents in the local and national surveys, respectively. Used to assign to and construct instruments for choice occasions by origin location.

2. [respondents_local.csv](/data/dwh_replication/respondents_local.csv) and [respondents_national.csv](/data/dwh_replication/respondents_national.csv): contain respondent demographic location information for all respondents in the local and national surveys, respectively. Used to assign to and construct instruments for choice occasions by origin location.

3. [trips_local.csv](/data/dwh_replication/trips_local.csv) and [trips_national.csv](/data/dwh_replication/trips_national.csv): contain information on shoreline recreation trips for respondents in the local and national surveys, respectively.  Used to assign to and construct instruments for choice alternatives by recreation site.

4. [trips_demos.txt](/data/dwh_replication/trips_demos.txt): contains the final, processed dataset used by English et al. (2018) alongside the travel cost dataset to estimate the shoreline recreation travel cost model. Includes chosen alternatives for each choice occasion as well as the necessary fields to create observation weights as discussed in the main text and in English et al. (2018).  I similarly use this alongside the travel cost dataset discussed below and the travel cost instruments discussed in the main text to implement the English et al. (2018) empirical exercise.

5. [tc1.txt](/data/dwh_replication/tc1.txt): contains the matrix of origin-destination travel cost estimates calculated by English et al. (2018). See English et al. (2018) for a discussion of the construction of these travel costs.

6. [local_caseids.txt](/data/dwh_replication/local_caseids.txt): contains the list of unique identifiers for all respondents from the local survey used in estimating the travel cost model using only this subsample (results reported in Appendix Tables A3 and A4).

#### Cost instrument data
I discuss cost instrument selection in Section 5.3 of the main text.  Travel cost instruments come from a variety of sources, including:

1. [acs5yr_2012.csv](/data/dwh_replication/acs5yr_2012.csv): contains a set of demographic variables aggregated to the Core Based Statistical Area (CBSA) level from the US Census Bureau's [2012 American Community Survey (ACS) 5-year Estimates](https://www.census.gov/data/developers/data-sets/acs-5year/2012.html).  These demographic variables include CBSA population, median household income, labor force size, number of employed individuals, and property tax rates.  These data are initially downloaded in their raw format directly from the US Census Bureau's API using the `tidycensus` package in R (Walker and Herman, 2024).

2. [annual_taxes_cleaned.dta](/data/dwh_replication/annual_taxes_cleaned.dta): contains raw time series data on state-level gasoline and other motor fuel excise taxes. Blake Barr provided excellent research assistance constructing measures of state motor fuels taxes.  NOTE: DO NOT REUSE WITHOUT PROPER ATTRIBUTION.  I request that you please contact me prior to reuse of these data to ensure the research assistance is properly attributed.

3. [eia_intl_crude_prod.txt](/data/dwh_replication/eia_intl_crude_prod.txt), [eia_refiner_costs.txt](/data/dwh_replication/eia_refiner_costs.txt), [oil_struc_shocks.csv](/data/dwh_replication/oil_struc_shocks.csv): contain data related to the construction of structural oil supply and demand shocks following Kilian (2009). [eia_intl_crude_prod.txt](/data/dwh_replication/eia_intl_crude_prod.txt) contains monthly global crude oil production data from 1973 to 2023 from the US Energy Information Administration (EIA) (EIA series id: "INTL.57-1-WORL-TBPD.M"). [eia_refiner_costs.txt](/data/dwh_replication/eia_refiner_costs.txt) contains monthly crude oil acquisition costs for US refiners from 1974 to 2023 obtained from the EIA (EIA series id: "PET.R1300____3.M"). To implement the method of Kilian (2009), I supplement these data with an index of global real economic activity maintained by the [Federal Reserve Bank of Dallas](https://www.dallasfed.org/research/igrea).  I obtain this index directly from the Federal Reserve Bank of St. Louis's API using the `fredr` package in R (Boysel and Vaughan, 2023). The file [oil_struc_shocks.csv](/data/dwh_replication/oil_struc_shocks.csv) contains the resulting estimated structural crude oil market shocks obtained following Kilian (2009).

4. [ZIP_CBSA_092012](/data/dwh_replication/ZIP_CBSA_092012.csv): contains a Zip code to CBSA crosswalk for 2012, which I obtain from the [US Department of Housing and Urban Development](https://www.huduser.gov/portal/datasets/usps_crosswalk.html).

5. [eia_wti.csv](/data/dwh_replication/eia_wti.csv), [nhgis0028_ds201_20135_zcta.csv](/data/dwh_replication/nhgis0028_ds201_20135_zcta.csv), and [respondents_commutes.csv](/data/dwh_replication/respondents_commutes.csv): all contain legacy cost instrument data used in earlier drafts of this paper and should therefore be viewed as deprecated.

#### Intermediate and final estimation data
These subfolders contain intermediate or final objects related to estimation of the recreation demand models as part of the English et al. (2018) empirical exercise.  In particular

1. [`/zmat/`](/data/dwh_replication/zmat/): contains wide format cost instrument matrices stored as CSV files for use in implementing the two-stage control function estimator as part of the English et al. (2018) empirical exercise. These are stored in wide format (i.e., a choice occasion-by-recreation site matrix) to align with the format of the final estimation objects used by English et al. (2018).  Note: this subfolder is empty in the final replication package to reduce storage requirements as these files are large in size.  Running the replication code as described below will populate this subfolder with the necessary wide format cost instrument matrices.

2. [`/est/`](/data/dwh_replication/est/): contains parameter estimates for all main and robustness specifications included in the English et al. (2018) empirical exercise.

3. [`/boot_est/`](/data/dwh_replication/boot_est/): contains parameter estimates for all bootstrap iterations for all main and robustness specifications included in the English et al. (2018) empirical exercise. These files are used to construct parameter covariances and standard errors.

## Computational Requirements

The code has been tested on the following machines:

- Dell Precision 5860 Tower with Intel ZXeon w5-2465X (3.10 GHz) processor with 16 cores, 192 GB of installed RAM, and 952 GB of storage running Windows 11.
- Apple MacBook Pro with Intel Core i5 (2.30 GHz) processor with 8 cores, 16 GB of isntalled RAM, and 500 GB of storage running  macOS Monterey 12.7.2.

I expect this package to run successfully on systems with similar specifications. Runtime will depend on the specifications of the system, particularly since much of the computationally-intensive steps written in Julia depend on multi-threading.  Note that replicators should adjust their enrivonment variable `JULIA_NUM_THREADS` to allow for multithreading in accordance with their system's specifications.

All code is written in R and Julia and has been tested with R version 4.4.1 and Julia version 1.10.4. The versions of the R packages are stored in the [renv](/renv.lock) lockfile. The versions of the required Julia dependencies are stored in the [Manifest.toml](/Manifest.toml) file.

## Instructions to Replicators

### 1. Installation of required software

Download and install version 4.4.1 of R and Julia 1.10.4.

### 2. Replicate the results of the paper

After installing the necessary software above, navigate to the project directory and exectute the [`run_all.jl`](run_all.jl) file.  This file will execute the code of the paper, starting with the code necessary to replicate the simulation results followed by the code necessary to replicate the empirical exercise.  This master file is written in Julia and uses the Julia package `RCall` to call all R scripts. The master file can be called using the following command line prompt in a Windows command terminal:

```
julia run_all.jl
```

You can also open a Julia REPL in your IDE of choice and execute the [`run_all.jl`](run_all.jl) script. Executing the [`run_all.jl`](run_all.jl) script will create all of the tables and figures of the paper. The results will be saved within the [`/output/`](/output/) folder, which contains separate subfolders for the [simulation results](/output/simulations/) and [empirical exercise results](/output/dwh_replication/).

The master file contains two commands which will activate and install the necessary R and Julia package dependencies:

```         
#   Activate project
using Pkg
Pkg.activate(".")

#   Load R packages:
R"renv::restore()"  ## Enter "y" when prompted
```

Follow all prompts to ensure proper installation of necessary R and Julia packages. 

## Source Code
All source code is stored within the [`/src/`](/src/) folder. The code is stored separately by language into subfolders [`/src/R/`](/src/R) and [`/src/julia/`](/src/julia/). Within each subfolder, files associated with the simulation results are prefixed by s\_ and files associated with the empirical results are prefixed by e\_.

The code will work if it is run in the order shown below; however, note that the above commands that will ensure all R and Julia dependencies are installed are only contained within the [`run_all.jl`](run_all.jl) file.  As a result, I encourage you to follow the automated replication process using the [`run_all.jl`](run_all.jl) script.

### Simulation Code

1. [`/julia/s_simulations.jl`](/src/julia/s_simulations.jl): contains all user-written functions necessary to simulate choice datasets and implement simulation results.

2. [`/julia/s_run_simulations.jl`](/src/julia/s_run_simulations.jl): contains commands to generate simulated datasets and implement all simulation results reported in the paper using the user-written functions contained within `s_simulations.jl` (loads these functions into the workspace). Saves intermediate data and estimation objects in [`/data/simulations/`](/data/simulations/). 

3. [`/R/s_simulation_figs.R`](/src/R/s_simulation_figs.R): contains code to generate all final tables and figures summarizing simulation results and stores these final outputs in [`/output/simulations/`](/output/simulations/). This script generates Figure 2 ([fig_dgp1.jpg](/output/simulations/fig_dgp1.jpg) and [fig_wtp1.jpg](/output/simulations/fig_wtp1.jpg)), Figure 3 ([fig_dgp2.jpg](/output/simulations/fig_dgp2.jpg) and [fig_wtp2.jpg](/output/simulations/fig_wtp2.jpg)), Figure 4 ([fig_wtp3.jpg](/output/simulations/fig_wtp3.jpg)), Figure A1 ([fig_rc1.jpg](/output/simulations/fig_rc1.jpg)), Figure A2 ([fig_rc2.jpg](/output/simulations/fig_rc2.jpg)), Figure A3 ([fig_rc3.jpg](/output/simulations/fig_rc3.jpg)), Figure A4 ([fig_rc4.jpg](/output/simulations/fig_rc4.jpg)), Figure A5 ([fig_rc5.jpg](/output/simulations/fig_rc5.jpg)), Table 1 ([tab_sim_est.tex](/output/simulations/tab_sim_est.tex)), and Table A1 ([tab_sim_est_groupfe.tex](/output/simulations/tab_sim_est_groupfe.tex)).

### Empirical Exercise Code
The code associated with the empirical exercise can be differentiated into two categories: data construction and recreation demand estimation code.

#### Data construction code:

1. [`/R/e_gen_oil_vars.R`](/src/R/e_gen_oil_vars.R): contains code to generate structural oil demand and supply shocks following Kilian (2009). Takes as input the [eia_intl_crude_prod.txt](/data/dwh_replication/eia_intl_crude_prod.txt), [eia_refiner_costs.txt](/data/dwh_replication/eia_refiner_costs.txt), and Dallas FED index of real economic activity and saves the resulting structural shocks in [oil_struc_shocks.csv](/data/dwh_replication/oil_struc_shocks.csv). Generates Figure A7 ([fig_shock_aggdem_wide.jpg](/output/dwh_replication/fig_shock_aggdem_wide.jpg), [fig_shock_oildem_wide.jpg](/output/dwh_replication/fig_shock_oildem_wide.jpg), and [fig_shock_oilsupply_wide.jpg](/output/dwh_replication/fig_shock_oilsupply_wide.jpg)).

2. [`/R/e_gen_gas_taxes.R`](/src/R/e_gen_gas_taxes.R): contains code to process source time series data on state-level motor vehicle excise taxes contained in [annual_taxes_cleaned.dta](/data/dwh_replication/annual_taxes_cleaned.dta) and stores the results in [gas_tax_delta.csv](/data/dwh_replication/gas_tax_delta.csv), which I omit from the replication package to conserve storage space (running this file will generate it). Generates Figure A6 ([fig_gas_tax.jpg](/output/dwh_replication/fig_gas_tax.jpg)).

3. [`/R/e_gen_dist_mat.R`](/src/R/e_gen_dist_mat.R): contains code to generate the origin-destination distance matrix used in constructing cost instruments. Uses as inputs [trips_local.csv](/data/dwh_replication/trips_local.csv), [respondents_local.csv](/data/dwh_replication/respondents_local.csv), [respondents_national.csv](/data/dwh_replication/respondents_national.csv), [origins_local.csv](/data/dwh_replication/origins_local.csv), and [origins_national.csv](/data/dwh_replication/origins_national.csv). Saves the constructed distance matrix as [distance_mat.csv](/data/dwh_replication/distance_mat.csv), which I omit from the replication package to conserve storage space (running this file will generate it).  Generates Figure 5 in the main text ([fig_dist_freq.jpg](/output/dwh_replication/fig_dist_freq.jpg) and [fig_site_map.jpg](/output/dwh_replication/fig_site_map.jpg)). Note that a basemap is downloaded using the R package `ggspatial` and the basemap tiles are stored in the folder [`/rosrm.cache/`](/rosm.cache/) (Dunnington, 2023).

4. [`/R/e_gen_census_vars.R`](/src/R/e_gen_census_vars.R): contains code to download the CBSA-level 5-year ACS estimates for demographic variables used as cost instruments via the `tidycensus` package (Walker and Herman, 2024). Saves the resulting data as [acs5yr_2012.csv](/data/dwh_replication/acs5yr_2012.csv).

5. [`/R/e_gen_commuting_time.R`](/src/R/e_gen_commuting_time.R): contains code to process data on Zip code level average commuting times.  This is a legacy file from earlier drafts as the resulting commuting time field is no longer used as a travel cost instrument.

6. [`/R/e_construct_final_inst.R`](/src/R/e_construct_final_inst.R): contains code to construct wide-format cost instrument matrices as well as code to estimate versions of first stage travel cost regressions reported in the paper. This code takes as input [distance_mat.csv](/data/dwh_replication/distance_mat.csv), [tc1.txt](/data/dwh_replication/tc1.txt), [trips_demos.txt](/data/dwh_replication/trips_demos.txt), [respondents_local.csv](/data/dwh_replication/respondents_local.csv), [respondents_national.csv](/data/dwh_replication/respondents_national.csv), [eia_wti.csv](/data/dwh_replication/eia_wti.csv) (deprecated), [oil_struc_shocks.csv](/data/dwh_replication/oil_struc_shocks.csv), [respondents_commutes.csv](/data/dwh_replication/respondents_commutes.csv) (deprecated), [ZIP_CBSA_092012](/data/dwh_replication/ZIP_CBSA_092012.csv), [gas_tax_delta.csv](/data/dwh_replication/gas_tax_delta.csv), and [acs5yr_2012.csv](/data/dwh_replication/acs5yr_2012.csv). Generates the wide-format cost instrument matrices contained in the subfolder [`/data/dwh_replication/zmat/`](/data/dwh_replication/zmat/). Also estimates the first-stage travel cost regression estimates reported in Table 2 (saved as [dwh_first_stage_main.tex](/output/dwh_replication/dwh_first_stage_main.tex)) and the first-stage travel cost regression estimates reported in Table A3 (saved as [dwh_first_stage_local.tex](/output/dwh_replication/dwh_first_stage_local.tex)). Note that I construct these tables of first-stage estimates here rather than in the Julia scripts that implement the two-stage control function (discussed below) as the first-stage estimation is nested within the complete two-stage procedure, which at present does not store the first stage results (but could with some minor tweaks to the code). The first-stage point estimates obtained here are identical to those found in the complete two-stage procedure.

#### Recreation demand estimation code:

1. [`/julia/e_nmnl.jl`](/src/julia/e_nmnl.jl): contains all user-written functions necessary to implement all estimators of the nested multinomial logit models of recreation demand outlined in the empirical exercise in Section 5 of the main text of the paper as well as the robustness results included in the appendix. Includes functions to estimate lost user day values as outlined in Appendix C. Each of the below scripts calls these functions.

2. [`/julia/e_est_nmnl_main.jl`](/src/julia/e_est_nmnl_main.jl): contains code that estimates the main empirical results reported in Table 3 within Section 5 of the main paper. Specifically, implements the baseline and two-stage control function estimators of the model outlined in Section 5. Takes as inputs [trips_demos.txt](/data/dwh_replication/trips_demos.txt), [tc1.txt](/data/dwh_replication/tc1.txt), and the contents of the subfolder [`/data/dwh_replication/zmat/`](/data/dwh_replication/zmat/). Saves point estimates and estimates from bootstrap iterations to [`/data/dwh_replication/est/`](/data/dwh_replication/est/) and [`/data/dwh_replication/boot_est/`](/data/dwh_replication/boot_est), respectively.  Generates Table 3 in the paper ([dwh_second_stage_main.tex](/output/dwh_replication/dwh_second_stage_main.tex)).

3. [`/julia/e_est_nmnl_main_rc.jl`](/src/julia/e_est_nmnl_main_rc.jl): contains code that estimates robustness checks of the main empirical results of Section 5 of the main paper. Specifically, implements the baseline and two-stage control function estimators of the model outlined in Section 5 using alternative combinations of first-stage cost instruments. Takes as inputs [trips_demos.txt](/data/dwh_replication/trips_demos.txt), [tc1.txt](/data/dwh_replication/tc1.txt), and the contents of the subfolder [`/data/dwh_replication/zmat/`](/data/dwh_replication/zmat/). Saves point estimates and estimates from bootstrap iterations to [`/data/dwh_replication/est/`](/data/dwh_replication/est/) and [`/data/dwh_replication/boot_est/`](/data/dwh_replication/boot_est), respectively.  Generates Table A2 in the paper ([dwh_second_stage_rc.tex](/output/dwh_replication/dwh_second_stage_rc.tex)).

4. [`/julia/e_est_nmnl_local.jl`](/src/julia/e_est_nmnl_local.jl): contains code that estimates robustness checks of the main empirical results of Section 5 of the main paper. Specifically, implements the baseline and two-stage control function estimators of the model outlined in Section 5 using just the local respondents. Takes as inputs [trips_demos.txt](/data/dwh_replication/trips_demos.txt), [tc1.txt](/data/dwh_replication/tc1.txt), [local_caseids.txt](/data/dwh_replication/local_caseids.txt), and the contents of the subfolder [`/data/dwh_replication/zmat/`](/data/dwh_replication/zmat/). Saves point estimates and estimates from bootstrap iterations to [`/data/dwh_replication/est/`](/data/dwh_replication/est/) and [`/data/dwh_replication/boot_est/`](/data/dwh_replication/boot_est), respectively.  Generates Table A4 in the paper ([dwh_second_stage_local.tex](/output/dwh_replication/dwh_second_stage_local.tex)).

## References

- Boysel S. and D. Vaughan. (2023). "fredr: An R Client for the 'FRED' API." R package version 2.1.0.9000, [https://github.com/sboysel/fredr](https://github.com/sboysel/fredr).
- Dunnington, D. (2023). "ggspatial: Spatial Data Framework for ggplot2." https://paleolimbot.github.io/ggspatial/, https://github.com/paleolimbot/ggspatial.
- English, E., R.H. von Haefen, J. Herriges, C. Leggett, F. Lupi, K. McConnell, M. Welsh, A. Domanski and N. Meade. (2018). "Estimating the value of lost recreation days from the Deepwater Horizon oil spill." Journal of Environmental Economics and Management, 91: 26–45.
- Kilian, L. (2009). "Not All Oil Price Shocks Are Alike: Disentangling Demand and Supply Shocks in the Crude Oil Market." American Economic Review, 99(3): 1053–1069.
- US Census Bureau. (2012). American Community Survey 2012. Retrieved from [https://data.census.gov](https://data.census.gov).
- US Department of Housing and Urban Development. (2024). HUD USPS ZIP Code Crosswalk Files. Retrieved from [https://www.huduser.gov/portal/datasets/usps_crosswalk.html](https://www.huduser.gov/portal/datasets/usps_crosswalk.html).
- Walker, K. and M. Herman. (2024). "tidycensus: Load US Census Boundary and Attribute data as 'tidyverse' and 'sf'-Ready Data Frames." R package version 1.6.6, [https://walker-data.com/tidycensus/](https://walker-data.com/tidycensus/)





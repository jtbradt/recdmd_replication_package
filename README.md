# Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models

## Author: Jacob Bradt

The code and data contained in this directory replicate 'Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models."

> **Abstract:** Conventional applications of recreation demand models likely suffer from two standard challenges with demand estimation, namely omitted variables bias and measurement error. Idiosyncratic prices in the form of individual-level travel costs can exacerbate these two challenges: the potential for non-random selection into travel costs through residential sorting and the difficulty of observing individual-level travel costs both work to bias traditional model estimates. I demonstrate the magnitude of this potential bias in conventional estimates of recreation demand models. I provide a relatively simple instrumental variables approach to address these two empirical challenges that substantially outperforms traditional estimates in numerical simulations. Replicating English et al. (2018), I find that accounting for potential selection into travel costs and measurement error through the instrumental variables approach decreases estimates of the welfare costs of the 2010 Deepwater Horizon oil spill by 12 percent.

The code in this replication package constructs the necessary simulated and empirical datasets used in the simulation and English et al (2018) replication exercises contained in the paper 'Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models.'  Subsequent code uses these datasets to estimate the main results of the paper, generating 11 of 12 figures and 7 of 7 tables contained within the main text and online appendices.  

All code is written in R and Julia. I have tested this replication package using both Windows and Mac operating systems.

Computation of the full set of simulations reported in the text is time-intensive.  As a result, these replication files reduce the number of Monte Carlo simulations from 1000 to 100 for each set of simulated results.  Similarly, standard errors reported in the English et al. (2018) replication exercise in the main text require a bootstrap procedure, which is computationally intensive.  I also reduce the number of bootstrap iterates from 200 to 20 in these replication results.  Users interested in obtaining the full results from the paper are welcome to adjust the number of simulations and bootstrap iterations.

With these reductions in computational burden, the replication package should take around 16 hours when run on a machine with 192GB of memory and 16 cores.

### Data Availability and Provenance Statements

-   [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

-   [x] All data **can be made** publicly available.

### Data Description

The data in this repository are organized into two separate folders, one each for the simulation results and the English et al. (2018) replication exercise, respectively.  The intermediate and final simulated data generated as part of the Monte Carlo analysis are stored in [data/simulations/](/data/simulations/).  The raw, intermediate, and final data used as part of the English et al. (2018) replication exercise are stored in [data/dwh_replication](/data/dwh_replication/).

#### [Simulation Data](/data/simula/)

1. [sim_data_sim1to6.csv](/data/simulations/sim_data_sim1to6.csv): contains single simulated data set for each of the data generating processes defined by Simulations 1-6 in the text.
2. [est_wtp_sim1.csv](/data/simulations/est_wtp_sim1.csv): contains estimated willingness-to-pay (WTP) statistics across all simulated datasets for the data generating process defined by Simulation 1 in the text.  Includes WTP estimates based on each of the three alternative estimators discussed in the text---the baseline multinomial logit estimator, the baseline multinomial logit estimator with origin group fixed effects, and two-stage control function estimator.  Similarly defined CSV files contain analogous WTP estimates for Simulations 2 through 6.
3. [est_wtp_rc1.csv](/data/simulations/est_wtp_rc1.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different degrees of within-recreator origin group correlation. The results from this robustness check are reported in Appendix Figure A1.
4. [est_wtp_rc2.csv](/data/simulations/est_wtp_rc2.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different recreator origin group sizes. The results from this robustness check are reported in Appendix Figure A2.
5. [est_wtp_rc3.csv](/data/simulations/est_wtp_rc3.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different forms of nonlinearity in the unobservable, idiosyncratic preference term. The results from this robustness check are reported in Appendix Figure A3.
6. [est_wtp_rc4.csv](/data/simulations/est_wtp_rc4.csv): contains estimated WTP statistics across all simulated datasets for the robustness check that implements the Simulation 2 data generating process with different forms of nonlinearity in the cost instrument. The results from this robustness check are reported in Appendix Figure A4.
7. [est_wtp_rc_weak.csv](/data/simulations/est_wtp_rc_weak.csv) and [est_wtp_rc_weak_fstat.csv](/data/simulations/est_wtp_rc_weak_fstat.csv): contains results from the robustness check that implements the Simulation 2 data generating process with different degrees of instrument relevance. The results from this robustness check are reported in Appendix Figure A5.
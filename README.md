# Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models

## Author: Jacob Bradt

The code in this replication package constructs the necessary simulated and empirical datasets used in the simulation and English et al (2018) replication exercises contained in the paper 'Hotelling Meets Wright: Spatial Sorting and Measurement Error in Recreation Demand Models.'  Subsequent code uses these datasets to estimate the main results of the paper, generating 11 of 12 figures and 7 of 7 tables contained within the main text and online appendices.  

All code is written in R and Julia. I have tested this replication package using both Windows and Mac operating systems.

Computation of the full set of simulations reported in the text is time-intensive.  As a result, these replication files reduce the number of Monte Carlo simulations from 1000 to 100 for each set of simulated results.  Similarly, standard errors reported in the English et al. (2018) replication exercise in the main text require a bootstrap procedure, which is computationally intensive.  I also reduce the number of bootstrap iterates from 200 to 20 in these replication results.  Users interested in obtaining the full results from the paper are welcome to adjust the number of simulations and bootstrap iterations.

With these reductions in computational burden, the replication package should take around 16 hours when run on a machine with 192GB of memory and 16 cores.

> **Abstract:** Conventional applications of recreation demand models likely suffer from two standard challenges with demand estimation, namely omitted variables bias and measurement error. Idiosyncratic prices in the form of individual-level travel costs can exacerbate these two challenges: the potential for non-random selection into travel costs through residential sorting and the difficulty of observing individual-level travel costs both work to bias traditional model estimates. I demonstrate the magnitude of this potential bias in conventional estimates of recreation demand models. I provide a relatively simple instrumental variables approach to address these two empirical challenges that substantially outperforms traditional estimates in numerical simulations. Replicating English et al. (2018), I find that accounting for potential selection into travel costs and measurement error through the instrumental variables approach decreases estimates of the welfare costs of the 2010 Deepwater Horizon oil spill by 12 percent.


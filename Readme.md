# Optimal energy management of offshore wind farms considering the combination of overplanting and dynamic rating
<img align="left" alt="Coding" width="200" src="https://www.cigre.org/userfiles/images/Events/2022/session-banner1.jpg">

  
This repository shares the MATLAB code and data for the conference paper 📋:\
Ildar Daminov, Anne Blavette, Salvy Bourguet, Didier Trichet, Guillaume Wasselynck, Laurent Dupont, Hamid Ben Ahmed, Thomas Soulard, and Pierre Warlop, "Optimal energy management of offshore wind farms considering the combination of overplanting and dynamic rating," in CIGRE session, Paris, France, 2022
  
  
## Paper's abstract
The electricity cost from offshore wind turbines has significantly declined over the past decade. However, in the context of energy transition, the further reduction in electricity generation costs remains a forefront subject. One possibility is to enable offshore wind farms to have a greater installed capacity than their transmission infrastructure, known as “overplanting”. This allows increasing their energy generation revenues while requiring some curtailments of their power output during the most energy-rich periods of the year. Also, export cables may have a high thermal inertia, so they can be fully loaded for several days before the cable reaches the maximum permissible temperature (usually 90℃ for XLPE cables). Hence, some Transmission System Operators (TSOs) have recently allowed offshore wind farm operators to export more power than the transmission infrastructure may deliver in the steady-state conditions. This is known as the «dynamic thermal rating» (DTR). DTR combined with overplanting can lead to a significant reduction in the Levelized Cost Of Energy (LCOE) of offshore wind farms. This paper examines the benefits of improving the electricity production commitment strategies against business-as-usual approaches based on the 50% quantile, called P50. In this perspective, a theoretical case is investigated, where the day-ahead and imbalance prices would be known in advance in order to define an optimized power production commitment. This theoretical case intends to define the upper bound on the annual revenue that could be gained by enhanced forecasts on both the production and the day-ahead and real-time energy prices. The study is based on a thermal model of a generic export cable, built per IEC standards 60287 and 60853-2, with data provided by the French TSO. The model was validated against simulation results provided by the French TSO. Historical day-ahead and imbalance prices are considered, as well as forecast and actual measured production profiles from offshore wind farms provided by the Belgian TSO. The source code developed in the project described in this paper, along with relevant data, and a related documentation will soon be provided in open-access repositories available from the project official webpage.
## How to run a code 
There are two ways how you may run this code: I. Do simulations yourself (but note that it takes few days) or II. Use precalculated data 
  
I. Launching all calculations yourself. This will reproduce the data in the paper but it would take several days:
1. Copy this repository to your computer 
2. Open the script main.m
3. Launch the script "main.m" by clicking on the button "Run" (usually located at the top of MATLAB window).\
As alternative, you may type ```main``` 
in Command Window to launch the entire script. 


II. Using the precalculated data to reproduce the particular figure : 
1. Copy this repository to your computer 
2. Open the script the Creating_Figures.m
3. Find the section (Figure XX) corresponding to the Figure you would like to reproduce. 
4. Put the cursor at any place of this section and click on the button "Run Section" (usually located at the top of MATLAB window)


## Files description

Principal scripts:
* main.m - the script which launches all calculations at a computer. Note that the entire calculates may take several days (due to neccesity to launch heavy and numerous optimization)
* Creating_Figures.m - this script repoduces the figures from the conference paper by using the precalculated data. 

Additional functions: 
* IEC60853_2_RTE_15min.m - a thermal model of submarine export cable 225 kV according to the standard IEC 60853-2 (15 min is a time resolution of input data)
* losses_partial_transients.m - this script estimates losses in a export cable. Also, it calculates partial transient (p 33 IEC 60853-2) of cable and its environment  
* losses_partial_transients_new.m - the same as losses_partial_transients.m but it requires 2 additional inputs (R_AC,Lambda1_new). These two parameters are changing as a function of temperature
* temperature_correction_on_losses_new.m -this function adjusts the final profile of temperature rises as a function of losses change due to a temperature variations.
* power_curtailement.m - this function estimates power curtailements of offshore wind farm if cable's limit of temperature (90 degC) is exceeded
* I_admissible.m - This function finds the admissible current in accordance with IEC 60287-1-1, p11
* I2P.m - this function converts I (in A) to power values (in MVA) for submarine cable 225 kV
* P2I.m - this function converts power S (in MVA) to current I (in A) for submarine cable 225 kV

Initial data:
* main_simulations_case1.mat - precalculated results for the case 1 (when Static Thermal Rating is used for cable constraints)
* main_simulations_case2.mat - precalculated results for the case 2 (when Dynamic Thermal Rating is used for cable constraints)
* data_DA_IMB_years_2015 2020.mat - day-ahead and imbalance prices for 2015-2020 (all at 15-min resolution). Source: [ENTSO-E. Transparency Platform](https://transparency.entsoe.eu/)
* Elia_Fev01_2012_Jan12_2017_preloads.mat - preload power profiles of offshore wind farm (2012-2017). Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* Elia_Fev01_2012_Jul28_2021_monitored_capacity_offshore.mat - monitored installed capacity of Belgian offshore wind farm (2012-2021). Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* Elia_Jan13_2016_Jan13_2017_powers.mat - Annual power profile of Belgian offshore wind farms. Both forecast and measurements. Source: [Elia](https://www.elia.be/en/grid-data/power-generation/wind-power-generation?csrt=6075160236430889381)
* ENTSO_E_Jan13_2016_Jan13_2017_prices.mat - Annual day-ahead and imbalance prices in France according to ENTSO-e. Source: [ENTSO-E. Transparency Platform](https://transparency.entsoe.eu/)
* theta_VECTOR_ISGT_main_cases_339MW_871A.mat - Temperature rises of export cable (after current 871 A was applied at preload intervals)



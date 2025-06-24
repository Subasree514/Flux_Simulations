## Codes to simulate the RS model integrated plant metabolic core model
Reactive Species (RS) model depicting various reactive species and antioxidant enzymatic reactions suitable for plant metabolism was constructed [not published yet]. RS model previously reconstructed for human metabolism was used as a template to reconstruct the plant RS model Sridhar, Subasree, et al.  _Metabolic Engineering_ (2023)

#### Codes to update the metabolic models with total demand reactions for the reactive species
https://github.com/Subasree514/Flux_Simulations/tree/trial/RS%20demands contains the codes to add demand reactions for
  1. Core metabolic model [not published yet] - core_demands.py
  2. Extended core metabolic model [not published yet] - Extended_core_demands.py

#### Codes to add generic and specific constraints to the models
https://github.com/Subasree514/Flux_Simulations/tree/trial/Constraints
  1. limiting_nutrient_analysis - Constraints added to simulate autotro[hic conditions
  2. exp_val_core.ipynb - Constraints added to antioxidant enzymes in core model
  3. exp_val.ipynb - Constraints added to antioxidant enzymes in extended core model

#### Codes to analyse the activity of reactions associated with specific metabolites at the equidistant points on the pareto plot
https://github.com/Subasree514/Flux_Simulations/tree/trial/Budget%20plots
  1. budget_H2O2.ipynb - Hydrogen peroxide associated reactions on the pareto plot
  2. budget_atp.ipynb - ATP producing/consuming reaction differences on the pareto curve
  3. budget_co2.ipynb - Carbondioxide associated reactions

biomass_rescaling.ipynb contains the code to rescale the biomass reaction to account for 1 g/mol per unit flux through it.

heatmap.ipynb contains the code to visually represent the variations in the fluxes through the reactions at the three chosen equidistant points on the pareto curve

rxns_from_compartment.py contains the code to find the number of reactions in each compartment

FVA_Copy_1 contains the code to get the flux spans of reactions along with reaction description

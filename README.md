## Codes to simulate the RS model integrated plant metabolic core model
Reactive Species (RS) model depicting various reactive species and antioxidant enzymatic reactions suitable for plant metabolism was constructed [not published yet]. RS model previously reconstructed for human metabolism was used as a template to reconstruct the plant RS model Sridhar, Subasree, et al.  _Metabolic Engineering_ (2023)

#### Codes to update the metabolic models with total demand reactions for the reactive species
https://github.com/Subasree514/Flux_Simulations/tree/trial/RS%20demands contains the codes to add demand reactions for
  1. Core metabolic model [not published yet] - core_demands.py
  2. Extended core metabolic model [not published yet] - Extended_core_demands.py

#### Codes to add generic and specific constraints to the models
https://github.com/Subasree514/Flux_Simulations/tree/trial/Constraints
  1. limiting_nutrient_analysis - Constraintsadded to simulate autotro[hic conditions
  2. exp_val_core.ipynb - Constraints added to antioxidant enzymes

#### Codes to analyse the activity of reactions associated with specific metabolites at the equidistant points on the pareto plot
https://github.com/Subasree514/Flux_Simulations/tree/trial/Budget%20plots
budget_H2O2.ipynb - Hydrogen peroxide associated reactions on the pareto plot
budget_atp.ipynb - ATP producing/consuming reaction differences on the pareto curve
budget_co2.ipynb - Carbondioxide associated reactions

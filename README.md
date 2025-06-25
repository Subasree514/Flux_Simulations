## Codes to simulate the RS model integrated plant metabolic core model
Reactive Species (RS) model depicting various reactive species and antioxidant enzymatic reactions suitable for plant metabolism was constructed [not published yet]. RS model previously reconstructed for human metabolism was used as a template to reconstruct the plant RS model Sridhar, Subasree, et al.  _Metabolic Engineering_ (2023)

#### Codes to update the metabolic models with total demand reactions for the reactive species
https://github.com/Subasree514/Flux_Simulations/tree/trial/RS%20demands contains the codes to add demand reactions for
  1. Core metabolic model [core_model_final.xml - not published yet] - core_demands.py. The model is now called as beta_day_DM.xml
  2. Extended core metabolic model [extended_core_model.xml - not published yet] - Extended_core_demands.py. The model is now called as beta_day_RS_DM_r.xml. extended_core_model.xml was constructed by merging the core_model_final.xml and the RS model for plants
     
#### Codes to add specific constraints to the models
Experimental values of activities of antioxidant enzymes are convereted into the bounds of the reespective reactions in the models. Experimental values for total demand demand reactions of non-enzymatic antioxidants like oxidised and reduced glutathione, ascorbate and dehydroascorbate were also added from the studies of (a)biotic stresses on plants
- https://github.com/Subasree514/Flux_Simulations/tree/trial/Constraints

  1. exp_val_core.ipynb - Constraints added to antioxidant enzymes in core model [beta_day_DM.xml]. The constrained model is now called as beta_antiox_core_dm.xml
  2. exp_val.ipynb - Constraints added to antioxidant enzymes in extended core model [beta_day_RS_DM_r.xml]. The constrained model is now called as beta_antiox_dm.xml

#### Codes to analyse the activity of reactions associated with specific metabolites at the equidistant points on the pareto plot in the constrained models
https://github.com/Subasree514/Flux_Simulations/tree/trial/Budget%20plots
  1. budget_H2O2.ipynb - Hydrogen peroxide associated reactions on the pareto plot
  2. budget_atp.ipynb - ATP producing/consuming reaction differences on the pareto curve
  3. budget_co2.ipynb - Carbondioxide associated reactions

biomass_rescaling.ipynb contains the code to rescale the biomass reaction to account for 1 g/mol per unit flux through it, to the models beta_day_RS_DM_r.xml and beta_day_DM.xml

heatmap.ipynb contains the code to visually represent the variations in the fluxes through the reactions at the three chosen equidistant points on the pareto curve

#### Archived codes
https://github.com/Subasree514/Flux_Simulations/blob/trial/RS%20demands/rename_core.py contains the codes to update the group names to the core model

https://github.com/Subasree514/Flux_Simulations/blob/trial/RS%20demands/rename_extendedcore.py contains the codes to update the group names to the extended core model

https://github.com/Subasree514/Flux_Simulations/blob/trial/Constraints/limiting_nutrient_analysis)limiting_nutrient_analysis - Constraints added to simulate autotrophic conditions

rxns_from_compartment.py contains the code to find the number of reactions in each compartment

FVA_Copy_1.py contains the code to get the flux spans of reactions along with reaction description

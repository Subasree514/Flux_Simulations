# Importing packages
from __future__ import print_function
import cobra
import os
from os.path import join
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cobra.util.solver import linear_reaction_coefficients
from IPython.display import Markdown, display
from cobra import flux_analysis

c3_model=cobra.io.load_matlab_model(join("alpha_day_RS_DM"))
solution = flux_analysis.pfba(c3_model)
solution_frame=solution.to_frame()
"""
Define search parameters

Select metabolites to search
"""

#Defining the list of metabolites to search
budget_metabolites = ['ATP_c','ATP_m','ATP_p','ATP_x']

#Defining the compartments within the model
compartment = ["_c", "_m", "_p", "_x", "_e", "_i", "_v", "_l", "tx", "ss"]

#Searching for the metabolites in the model
#for met in c3_model.metabolites.query("ATP"):
    #if met.id[:3] == "[M]":
 #   budget_metabolites.append(met)

print(budget_metabolites)
"""
Define consuming and producing reactions
"""

#Defining list of reactions producing and consuming the metabolite
consumers = []
producers = []

#Add reactions to respective list and exclude transport reactions
for met in budget_metabolites:
    for reaction in c3_model.reactions:
        if met in reaction.reactants and reaction.id[-2:] in compartment:
            consumers.append(reaction.id)
        elif met in reaction.products and reaction.id[-2:] in compartment:
            producers.append(reaction.id)
"""
Correct consumption/production with regards to directionality
"""
def budget_plot(solution_frame, producers, consumers):
    
    #Get flux values from the simulation for metabolite consuming/producing reactions
    producers_df = solution_frame.loc[producers,:]
    consumers_df = solution_frame.loc[consumers,:]

    #Get values with negative flows: producing reactions with negative flow are consuming and vice-versa
    negative_producers = list(producers_df[producers_df["fluxes"] < 0].index)
    negative_consumers = list(consumers_df[consumers_df["fluxes"] < 0].index)


    #Add reactions to correct list
    consumers.extend(negative_producers)
    producers.extend(negative_consumers)

    #Remove reactions with negative flux from old list
    def remove_items(test_list, item):
        res = [i for i in test_list if i != item]
        return res

    for item in negative_producers:
        producers = remove_items(producers, item)

    for item in negative_consumers:
        consumers = remove_items(consumers, item)

    """
    Get producing reactions and fluxes
    """

    #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
    producers_df = solution_frame.loc[producers,:]
    #Make all values positive (disregard directionality)
    print(producers_df)
    producers_df["fluxes"] = producers_df["fluxes"].abs()
    #Remove reactions with zero flux
    producers_df = producers_df[(producers_df.T != 0).all()]
    producers_df


    """
    Get consuming reactions and fluxes
    """

    #Get flux values from the simulation for metabolite consuming/producing reactions (correct list)
    consumers_df = solution_frame.loc[consumers,:]
    print(consumers_df)
    #Make all values positive (disregard directionality)
    consumers_df["fluxes"]  = consumers_df["fluxes"].abs()
    #Remove reactions with zero flux3_P_SERINE_
    consumers_df = consumers_df[(consumers_df.T != 0).all()]
    consumers_df

    """
    Concatenate producer and consumer dataframes
    """

    producers_df["Status"] = "Producer"
    consumers_df["Status"] = "Consumer"

    frame = [producers_df, consumers_df]

    all_reactions = pd.concat(frame)

    all_reactions["label"] = all_reactions.index

    #Export full dataframe to csv
    #all_reactions.to_csv("budget_plot.csv", header=True)

    #Correct Budget stoichiometry
    atp_stoi=[]
    for i in all_reactions['label']:
        for met in budget_metabolites:
            try:
                rxn=c3_model.reactions.get_by_id(i).get_coefficient(met)
                atp_stoi.append(abs(rxn))
            except KeyError:
                continue

    new_flux = all_reactions["fluxes"] * atp_stoi

    all_reactions.insert(3,'coefficient', atp_stoi)
    all_reactions.insert(4,'new_flux', new_flux)

    """
    Fluxes in consumption and production should be equal
    """

    #Sum the flux values
    print("Sum of fluxes: {}".format(all_reactions.groupby(["Status"]).new_flux.sum()))

    """
    Pick the colors - using random
    """

    import random
    import matplotlib.pyplot as plt


    #Defining the nº of colors
    number_of_colors = len(all_reactions.index)

    #Getting a list of colors
    random.seed(177)
    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                for i in range(number_of_colors)]

    #Getting list of reactions
    reaction_list = list(all_reactions.index)

    #Build color dictionary
    color_dict = {}
    for i in range(len(reaction_list)):
        color_dict[reaction_list[i]] = color[i]

    """
    Plot the pivot table and barplot
    """

    plt.style.use('fast')

    chart = all_reactions.pivot_table(index="Status", columns="label", values="new_flux")
    chart.plot.bar(rot = 0, stacked = True, legend = True, ylabel = "Flux", color = color_dict)
    plt.legend(loc='best', bbox_to_anchor=(1.0, 0.5, 0.5, 0.5), ncol = 2)
    #plt.title("Platoquinone Turnover in the Bundle Sheath Cell")
    figsize = [11, 11] #To prevent the cropping of the image
    #plt.savefig(plot_name)
    chart
budget_plot(solution_frame, producers, consumers)

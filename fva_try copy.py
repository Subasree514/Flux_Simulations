from __future__ import division, print_function, absolute_import
import csv
import pandas
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import os
import xml.etree.ElementTree as etree
import cobra
import numpy as np
from itertools import chain
from cobra.util import solver as sutil
from cobra.core.solution import get_solution
from optlang.symbolics import add, Zero
import pandas as pd
import os
from os.path import join
import matplotlib.pyplot as plt
from cobra.medium import minimal_medium
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
from cobra.flux_analysis import production_envelope
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model

model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
model_rs.solver='glpk'
print(model_rs.groups)
## Pareto function
# Pareto
objective1={''}
objective2={''}
pareto_range = (0.0, 1.001)  # for some reason you need to pick a number higher than 1).
pareto_step_size = 0.01
analysis_type = 'pareto'
metric = 'manhattan'
rxn2avoid = {''}
solver='gurobi'
constants = {'deltaC_CO2': 0.0055, 'D_H2O_0': 2.13E-05, 'D_CO2_0': 1.33E-05, 'mid_day': 6, 'deltaT': 2,
             'FeasTol': 1e-03, 'OptTol': 1e-03}
def pareto_analysis(model, objective1=objective1, objective2=objective2, pareto_range=pareto_range, metric=metric):
    reaction_obj1 = model.reactions.get_by_id(objective1)
    reaction_obj2 = model.reactions.get_by_id(objective2)
    result_list = []
    model.objective = {}
    reaction_obj1.objective_coefficient = 1
    solution = model.optimize()
    #print("\nSolving model (FBA) for determining objective 1 flux...")
    max_obj1 = dict(solution.fluxes)[objective1]
    #print("Max {0}: {1}".format(objective1, max_obj1))
    # change objective
    reaction_obj1.objective_coefficient = 0
    reaction_obj2.objective_coefficient = 1
    #print("\nSolving all iterations for Pareto frontier (FBA)...")
    for pareto in np.arange(pareto_range[0], pareto_range[1], pareto_step_size):
        if pareto == 1:
            reaction_obj1.lower_bound = max_obj1 * pareto  # * 0.999 # we need to add a bit of slack as the quadratic optimization is less accurate than the linear couterpart
        else:
            reaction_obj1.lower_bound = max_obj1 * pareto  # * 0.9999
        sol = model.optimize(objective_sense='maximize')
        # fix this minimal water loss value
        reaction_obj2.bounds = (sol.get_primal_by_id(objective2), sol.get_primal_by_id(objective2))
        if metric == 'manhattan':
            solution = cobra.flux_analysis.pfba(model)
            #solution = model.optimize()
            # print({'proline sink': solution['SK_PRO_c_06'], 'biomass 05': solution['Leaf_biomass_tx_05'], 'biomass 06': solution['Leaf_biomass_tx_06']})
            # solution.fluxes.to_excel(f'pareto_no_{pareto}.xlsx')
            result_list.append([pareto, solution[objective1], solution[objective2]])
            reaction_obj2.bounds = (0, 1000.0)
        elif metric == 'euclidean':

            # make copy because that is easier that reverting all the solver settings
            copy_model = model.copy()
            model.solver = solver

            FeasTol = float(constants['FeasTol'])
            OptTol = float(constants['OptTol'])

            copy_model.solver.configuration.tolerances.feasibility = FeasTol
            copy_model.solver.configuration.tolerances.optimality = OptTol

            rxnlist = [r for r in copy_model.reactions if r.id not in rxn2avoid]

            obj_vars = chain.from_iterable([r.flux_expression ** 2] for r in rxnlist)
            copy_model.objective = copy_model.problem.Objective(add(obj_vars), direction='min')

            print('\nSolving quadratic minimisation of sum of fluxes')
            #print(solver)
            solution = copy_model.optimize(objective_sense=None)
            result_list.append([pareto, solution[objective1], solution[objective2]])
        reaction_obj2.bounds = (0, 1000.0)
    return result_list
## Plots
#sol = core_model.optimize()
#print(core_model.summary(sol))
#print(core_model.reactions)
f1 = plt.figure(1)
## plot pareto plots
objective1 =  'H2O2_m_demand'
objective2 =  'RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p'
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
core_model=model_rs
#print(core_model.reactions.query('H2O2'))
rubisco = core_model.problem.Constraint(1 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
# Adding to model
core_model.add_cons_vars([rubisco])
##
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
data=pd.DataFrame(result_list)
plt.plot(data[1],data[2])
plt.xlabel('Nitric oxide demand')
plt.ylabel('Rubisco carboxylase')
plt.title("Nitric oxide vs. Rubisco carboxylase at Vc/Vo=1:1")
plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/no_co2_1_rs.pdf')
plt.show()
#with core_model:
#   core_model.reactions.get_by_id('DM_no[cell]').bounds = (0, 0)
#   sol = core_model.optimize()
#    print(core_model.summary(sol))
#   #result_list_0_rs=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#with core_model:
#   core_model.reactions.get_by_id('DM_HS_cell[cell]').bounds = (18,18)
#   sol = core_model.optimize()
#   print(core_model.summary(sol))
#   #result_list_half_rs=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#with core_model:
#   core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (150, 150)
#   sol = core_model.optimize()
#    print(core_model.summary(sol))
#    result_list_max_rs=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
#data=pd.DataFrame(result_list_max_rs)
#rint(data[2][0])
#print(data[2][100])
#plt.plot(data[1],data[2])
#plt.xlabel('Hydrogen peroxide demand')
#plt.ylabel('Biomass')
#plt.title("Hydrogen peroxide vs. Biomass reaction")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/h2o2_Biomass.pdf')
#plt.show() 
## important code (yet to analyse)
#maxs=[]
#mins=[]
#rxn_list=['FUMHYDR_RXN_c',	'ATP_pc',	'GLYCINE_AMINOTRANSFERASE_RXN_x',	'ALANINE_GLYOXYLATE_AMINOTRANSFERASE_RXN_x',	'MAL_CIT_vc',	'MAL_CIT_rev_vc',	'H_pc',	'Pi_PROTON_mc',	'H_mc',	'H_im',	'PYRUVATE_pc',	'F16ALDOLASE_RXN_p',	'NITRATE_vc',	'GLT_MAL_pc',	'2KG_MAL_pc',	'2.7.1.90_RXN_c',	'PGLUCISOM_RXN_c',	'F16ALDOLASE_RXN_c',	'G6P_Pi_pc',	'PGLUCISOM_RXN_p',	'PRO_mc',	'PRO_GLU_mc',	'CO2_mc',	'ISOCITDEH_RXN_c',	'NH3_mc',	'CARBAMATE_KINASE_RXN_p',	'MAL_xc',	'TRANSALDOL_RXN_p',	'MALATE_DEH_RXN_x',	'GLYCOLLATE_pc',	'MALATE_DEH_RXN_m',	'GLUTAMATE_SYNTHASE_NADH_RXN_p',	'GLUTAMINESYN_RXN_m',	'GLUTAMINESYN_RXN_p',	'unlProtHYPO_c',	'GLUTAMINESYN_RXN_c',	'GLN_GLU_mc',	'PEPDEPHOS_RXN_p',	'FRU_vc',	'GDPKIN_RXN_c',	'GLYCOGENSYN_RXN_p',	'PEPDEPHOS_RXN_c',	'GLY3KIN_RXN_p',	'SUCROSE_PHOSPHATASE_RXN_c',	'GLC_vc',	'RXN_1826_p',	'CIT_PROTON_vc',	'NADH_KINASE_RXN_c',	'SUCROSE_SYNTHASE_RXN_c',	'PSERTRANSAM_RXN_p',	'1.2.1.9_RXN_c',	'GLUTAMATE_DEHYDROGENASE_RXN_m',	'SEDOHEPTULOSE_BISPHOSPHATASE_RXN_p',	'PEPCARBOX_RXN_c',	'MAL_PROTON_vc',	'Glycerate_xc',	'HYDROXYPYRUVATE_REDUCTASE_RXN_NAD_x',	'RXN_13202_p',	'RXN_12486_c',	'2.7.7.13_RXN_c',	'ISOCITRATE_DEHYDROGENASE_NAD_RXN_m',	'ATPase_tx',	'INORGPYROPHOSPHAT_RXN_c',	'PEPCARBOXYKIN_RXN_c',	'F16BDEPHOS_RXN_p',	'GLUC1PADENYLTRANS_RXN_p',	'PROTONATP_rev_vc',	'6PFRUCTPHOS_RXN_c',	'2.7.7.34_RXN_c',	'GLC_PROTON_rev_vc',	'MANNPGUANYLTRANGDP_RXN_c',	'RXN_7703_c',	'6PFRUCTPHOS_RXN_p',	'RXN0_5224_c',	'MAL_PROTON_rev_vc',	'GLC_pc',	'PROTON_ATPase_c',	'F16BDEPHOS_RXN_c',	'SEDOBISALDOL_RXN_p',	'FORMATETHFLIG_RXN_p',	'NITRATE_PROTON_rev_vc',	'PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c',	'H_tx',	'PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p',	'RXN0_5114_p',	'FRU_PROTON_rev_vc',	'SUCROSE_PHOSPHATE_SYNTHASE_RXN_c',	'PROTON_PPi_rev_vc',	'PGLYCDEHYDROG_RXN_p',	'FORMYLTHFDEFORMYL_RXN_p',	'SERINE_GLYOXYLATE_AMINOTRANSFERASE_RXN_x',	'GLYCERATE_GLYCOLLATE_pc',	'FRUCTOKINASE_RXN_c',	'INORGPYROPHOSPHAT_RXN_p',	'CO2_pc',	'Mehler_Reaction_p',	'Photon_tx',	'RS_Plant_21_C',	'Photon_ep',	'RXN490_3650_p',	'1.2.1.13_RXN_p',	'RS_HOT5_CP',	'PHOSPHOGLUCMUT_RXN_c',	'ISOCITDEH_RXN_m',	'ORNITHINE_CITRULLINE_pc',	'L_CITRULLINE_pc',	'H2O_tx',	'H2O_ec',	'ASPAMINOTRANS_RXN_p',	'RS_Plant_2_P',	'RS_Plant_3_CP',	'H2O_xc',	'1.18.1.2_RXN_p',	'GLUTAMATE_SYNTHASE_FERREDOXIN_RXN_p',	'PYRROLINECARBDEHYDROG_RXN_NADP_m',	'THR_PROTON_vc',	'MALIC_NADP_RXN_p',	'HEXOKINASE_RXN_MANNOSE_c',	'GLN_PROTON_vc',	'MANNITOL_1_PHOSPHATASE_RXN_c',	'GLUTKIN_RXN_c',	'LEU_PROTON_rev_vc',	'L_ALPHA_ALANINE_PROTON_vc',	'L_ASPARTATE_PROTON_vc',	'L_ALPHA_ALANINE_PROTON_rev_vc',	'SUCROSE_PROTON_rev_vc',	'MET_PROTON_vc',	'PHE_PROTON_vc',	'HIS_PROTON_vc',	'1.1.1.255_RXN_c',	'1.1.1.39_RXN_m',	'CELLULOSE_SYNTHASE_UDP_FORMING_RXN_c',	'ARG_PROTON_vc',	'ASNSYNA_RXN_c',	'MANNOSE_6_PHOSPHATE_6_REDUCTASE_RXN_c',	'CONIFERIN_BETA_GLUCOSIDASE_RXN_c',	'L_ASPARTATE_PROTON_rev_vc',	'TRP_PROTON_rev_vc',	'THR_PROTON_rev_vc',	'ARG_PROTON_rev_vc',	'TRP_PROTON_vc',	'GLT_PROTON_vc',	'GLN_PROTON_rev_vc',	'HIS_PROTON_rev_vc',	'LYS_PROTON_rev_vc',	'CIT_PROTON_rev_vc',	'GLY_PROTON_vc',	'4_AMINO_BUTYRATE_PROTON_vc',	'UDPKIN_RXN_c',	'ILE_PROTON_vc',	'RXN_10773_c',	'GLUC1PURIDYLTRANS_RXN_c',	'CELLULOSE_SYNTHASE_GDP_FORMING_RXN_c',	'CYS_PROTON_vc',	'PRO_PROTON_vc',	'MALTODEXGLUCOSID_RXN_p',	'LYS_PROTON_vc',	'GLY_PROTON_rev_vc',	'SER_PROTON_vc',	'ASN_PROTON_vc',	'MET_PROTON_rev_vc',	'GLT_PROTON_rev_vc',	'VAL_PROTON_vc',	'SUCROSE_PROTON_vc',	'TYR_PROTON_vc',	'RS_Plant_28_C',	'RS_Plant_OXP1_C',	'PRO_PROTON_rev_vc',	'ASPARAGHYD_RXN_c',	'VAL_PROTON_rev_vc',	'SER_PROTON_rev_vc',	'MALIC_NADP_RXN_c',	'GLUCOKIN_RXN_c',	'CYS_PROTON_rev_vc',	'2.4.1.111_RXN_c',	'ASN_PROTON_rev_vc',	'LEU_PROTON_vc',	'PHE_PROTON_rev_vc',	'4_AMINO_BUTYRATE_PROTON_rev_vc',	'GLUCOKIN_RXN_p',	'TYR_PROTON_rev_vc',	'ILE_PROTON_rev_vc',	'Pi_xc',	'H_xc',	'O2_pc',	'PLASTOQUINOL_PLASTOCYANIN_REDUCTASE_RXN_p',	'O2_tx',	'O2_ec',	'ho2_rad_p_demand',	'DM_oh_rad[cell]',	'H2O2_p_demand',	'Ferredoxin_Plastoquinone_Reductase_p',	'oh_rad_p_demand',	'DM_HYDROGEN_PEROXIDE_cell[cell]',	'NADPH_Dehydrogenase_p',	'N_ACETYLGLUTPREDUCT_RXN_p',	'GLUTAMATE_N_ACETYLTRANSFERASE_RXN_p',	'ACETYLORNTRANSAM_RXN_p',	'ACETYLGLUTKIN_RXN_p',	'ORNITHINE_GLU_AMINOTRANSFORASE_RXN_m',	'SPONTPRO_RXN_m',	'GLY_xc',	'PHOSPHORIBULOKINASE_RXN_p',	'O2S_c_demand',	'DM_SUPER_OXIDE_cell[cell]',	'O2S_p_demand',	'RS_76_CP',	'MDA_Fd_Ascorbate_p',	'DM_oh1[cell]',	'oh_p_demand',	'DSERDEAM_RXN_c',	'6PGLUCONOLACT_RXN_p',	'6PGLUCONDEHYDROG_RXN_p',	'RXN0_5184_c',	'GLU6PDEHYDROG_RXN_p',	'MALTODEG_RXN_c',	'GLU6PDEHYDROG_RXN_c',	'RIBULP3EPIM_RXN_c',	'3.2.1.48_RXN_c',	'RXN_1841_v',	'RXN_14351_pc',	'RXN_1781_v',	'RXN_1827_p',	'6PGLUCONDEHYDROG_RXN_c',	'3.2.1.48_RXN_v',	'6PGLUCONOLACT_RXN_c',	'ACETATE_COA_LIGASE_RXN_p',	'RS_128_CP',	'RS_Plant_NAGS_CP',	'ACET_pc',	'Plastidial_ATP_Synthase_p',	'ARG_pc',	'ORNCARBAMTRANSFER_RXN_p',	'UREA_mc',	'RXN_2141_p',	'UREASE_RXN_c',	'ARGSUCCINSYN_RXN_p',	'ARGINASE_RXN_m',	'FUM_pc',	'ARGSUCCINLYA_RXN_p',	'THRESYN_RXN_p',	'ASPARTATEKIN_RXN_p',	'AMP_ATP_xc',	'HOMOSERKIN_RXN_p',	'INORGPYROPHOSPHAT_RXN_x',	'ASPARTATE_SEMIALDEHYDE_DEHYDROGENASE_RXN_p',	'MALSYN_RXN_x',	'THREONINE_ALDOLASE_RXN_c',	'ACETATE_COA_LIGASE_RXN_x',	'HOMOSERDEHYDROG_RXN_NAD_p',	'RXN66_3_c',	'THR_pc',	'RXN66_3_m',	'ACET_mc',	'RXN0_5330_NAD_mi',	'Plastoquinol_Oxidase_p',	'1.10.2.2_RXN_mi',	'SPONTPRO_RXN_c',	'PYRROLINECARBREDUCT_RXN_NADP_c',	'NADPHoxc_tx',	'RXN0_5330_NADP_mi',	'RS_27_CP',	'RXN0_5330_NAD_c',	'RS_74_CP',	'RS_46_CP',	'RXN0_5330_NADP_mc',	'RXN_14903_mi',	'H2O2_v_demand',	'1.10.2.2_RXN_mc',	'H2O2_c_demand',	'GLYC3PDEHYDROGBIOSYN_RXN_c',	'SUPEROX_DISMUT_RXN_p',	'RS_35_CP',	'NADPHoxp_tx',	'NADH_DEHYDROG_A_RXN_mc',	'RXN_12541_c',	'PSII_RXN_p',	'RS_Plant_H2O2_tr_V',	'RS_25_CP',	'NADH_DEHYDROG_A_RXN_mi',	'FERRIC_CHELATE_REDUCTASE_RXN_c',	'RXN_6883_mc',	'RXN_14903_m',	'Mitochondrial_ATP_Synthase_m',	'RXN0_5330_NADP_c',	'RXN_6883_mi',	'GLUTSEMIALDEHYDROG_RXN_c',	'NADPHoxm_tx',	'RXN0_5330_NAD_mc',	'SUPEROX_DISMUT_RXN_c',	'RXN0_5260_m',	'RS_46_C',	'OOH_c_demand',	'DM_OOH-[cell]',	'METHENYLTHFCYCLOHYDRO_RXN_c',	'1.2.1.2_RXN_m',	'METHYLENETHFDEHYDROG_NADP_RXN_c',	'FORMATE_pc',	'1.2.1.2_RXN_p',	'HOMOSERDEHYDROG_RXN_NADP_p',	'CITSYN_RXN_m',	'SUC_xc',	'SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mi',	'PYRUVDEH_RXN_m',	'ISOCIT_CLEAV_RXN_x',	'GLUTATHIONE_REDUCT_NADPH_RXN_p',	'RXN_3523_p',	'RXN_3522_p',	'1.8.5.1_RXN_p',	'RXN_9623_p',	'PHOSPHOLIPASE_A2_RXN_c',	'2.3.1.23_RXN_p',	'SERINE_O_ACETTRAN_RXN_c',	'LCYSDESULF_RXN_c',	'ACSERLY_RXN_c',	'ATP_CITRATE_PRO_S_LYASE_RXN_c',	'CITSYN_RXN_x',	'CIT_xc',	'RXN_3521_p',	'THIOREDOXIN_REDUCT_NADPH_RXN_p',	'L_ASCORBATE_PEROXIDASE_RXN_p',	'RS_Plant_PER1_CP',	'GLUTATHIONE_PEROXIDASE_RXN_p',	'SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mc',	'SUCCCOASYN_RXN_m',	'GCVMULTI_RXN_m',	'RIB5PISOM_RXN_p',	'SUCCINATE_COA_LIGASE_GDP_FORMING_RXN_m',	'1TRANSKETO_RXN_p',	'2OXOGLUTARATEDEH_RXN_m',	'PYRUVDEH_RXN_p',	'SO4_ec',	'1.8.4.9_RXN_p',	'SO3_p_demand',	'DM_SO3_cell[cell]',	'SULFATE_ADENYLYLTRANS_RXN_p',	'SO4_tx',	'SUCCINYL_COA_HYDROLASE_RXN_m',	'GLUTDECARBOX_RXN_c',	'GABATRANSAM_RXN_m',	'SUCCINATE_SEMIALDEHYDE_DEHYDROGENASE_RXN_m',	'ALANINE_AMINOTRANSFERASE_RXN_m',	'RXN_6902_m',	'RXN_6884_m',	'CYTOCHROME_C_OXIDASE_RXN_mi',	'CYTOCHROME_C_OXIDASE_RXN_mc',	'RXN66_1_c',	'RXN_7978_p',	'RS_Plant_DHAR_C',	'RS_Plant_GPX2_C',	'RXN_7985_p',	'RS_Plant_GR_C',	'RXN_7984_p',	'RS_Plant_APX_C',	'RXN0_1483_c',	'RXN_7979_p',	'MALONYL_COA_ACP_TRANSACYL_RXN_p',	'H2O2_x_demand',	'RXN0_5224_p',	'ACETYL_COA_CARBOXYLTRANSFER_RXN_p',	'RS_Plant_MDHAR_C',	'3_HYDROXBUTYRYL_COA_DEHYDRATASE_RXN_x',	'CO2_xc',	'DIHYDROPICRED_RXN_NAD_p',	'BHBDCLOS_RXN_x',	'DIAMINOPIMEPIM_RXN_p',	'DIHYDRODIPICSYN_RXN_p',	'ALLYSINE_DEHYDROG_RXN_c',	'LYS_pc',	'GLUTACONYL_COA_DECARBOXYLASE_RXN_x',	'1.5.1.9_RXN_c',	'GLUTACONYL_COA_mc',	'GLUTARYL_COA_DEHYDROG_RXN_m',	'1.5.1.8_RXN_c',	'RXN_7737_p',	'DIAMINOPIMDECARB_RXN_p',	'ACETYL_COA_ACETYLTRANSFER_RXN_x',	'2_KETO_ADIPATE_DEHYDROG_RXN_m',	'2_AMINOADIPATE_AMINOTRANSFERASE_RXN_c',	'DIHYDROPICRED_RXN_NADP_p',	'RXN_969_x',	'GPH_RXN_p',	'SULFITE_REDUCTASE_FERREDOXIN_RXN_p',	'HS_pc',	'RXN_13161_m',	'H2O2_m_demand',	'MercaptoPyruvateSulfurtransferase_m',	'CYSTEINE_AMINOTRANSFERASE_RXN_m',	'SULFITE_OXIDASE_RXN_m',	'RXN_7676_p',	'RXN_7678_NAD_p',	'RXN_7677_p',	'H2S_c_demand',	'RXN_7674_p',	'H2S_p_demand',	'DM_HS_cell[cell]',	'RXN_7678_NADP_p',	'RXN_7679_p',	'CATAL_RXN_x',	'SO3_m_demand',	'GLUTATHIONE_REDUCT_NADPH_RXN_m',	'L_ASCORBATE_PEROXIDASE_RXN_m',	'1.8.5.1_RXN_m',	'RS_Plant_GPX_M3',	'RXN_9958_NAD_m',	'RXN_6384_x',	'RXN_11832_p',	'3_HYDROXYPROPIONATE_DEHYDROGENASE_RXN_m',	'3_HYDROXY_PROPIONATE_xc',	'PROPCOASYN_RXN_x',	'HMBPP_synthesis_p',	'ISPH2_RXN_p',	'RS_Plant_CAT_M',	'RS_Plant_2_M',	'RS_56_M',	'RS_27_M',	'RS_Plant_1_M',	'oh_m_demand',	'RS_76_M',	'RS_25_M',	'RXN_9531_p',	'RXN_9532_p',	'PALMITATE_pc',	'RXN_9660_p',	'RXN_9528_p',	'RXN_9653_p',	'RXN_9533_p',	'RXN_9648_p',	'RXN_9659_p',	'RXN_9520_p',	'RXN_9650_p',	'RXN_9535_p',	'RXN_9658_p',	'4.2.1.61_RXN_p',	'RXN_9537_p',	'RXN_9661_p',	'Beta_Oxidation_x',	'RXN_9524_p',	'RXN_9662_p',	'RXN_9539_p',	'2.3.1.180_RXN_p',	'RXN_9652_p',	'RXN_9651_p',	'RXN_9527_p',	'4.2.1.58_RXN_p',	'RXN_9523_p',	'RXN_9655_p',	'RXN_9540_p',	'RXN_9654_p',	'RXN_9657_p',	'RXN_9536_p',	'RXN_9518_p',	'RXN_9514_p',	'RXN_9549_p',	'4.2.1.59_RXN_p',	'RXN_9663_p',	'RXN_9516_p',	'NO_p_demand',	'NITRITE_pc',	'FERREDOXIN_NITRITE_REDUCTASE_RXN_p',	'NO_m_demand',	'RS_Tr_no2',	'Nitrate_tx',	'Nitrate_ec',	'NITRATE_REDUCTASE_NADH_RXN_c',	'DM_no[cell]',	'Phytol_degradation_p',	'PROPIONYL_COA_mc',	'METHYLACYLYLCOA_HYDROXY_RXN_m',	'Phytol_biosynthesis_p',	'A_B_oxidation_x',	'ISOBUTYRYL_COA_xc',	'RXN0_884_p',	'3_HYDROXYISOBUTYRYL_COA_HYDROLASE_RXN_m',	'GGPP_biosynthesis_p',	'3_HYDROXYISOBUTYRATE_DEHYDROGENASE_RXN_m',	'CPD_14927_pc',	'MEPROPCOA_FAD_RXN_m',	'ISOCITDEH_RXN_x',	'RXN_11213_m',	'RS_2_M',	'RS_2_CP']
#for i in rxn_list[:10]:
#   result_list=pareto_analysis(core_model, objective1 = objective1, objective2=i, pareto_range = pareto_range, metric = metric)
#   data=pd.DataFrame(result_list)
#   if data[2][0] == data[2][100]:
#      maxs.insert(0,i)
       #maxs.extend([i])
        #maxs.append([i])
#print(maxs)
#print(type(maxs))



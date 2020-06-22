#Written by Keiran Rowell (z3374843) for use with PIEDA_mat.py for analysis of FMO data.
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
interactions_list = []
energies_list = sys.argv[1]
with open(energies_list) as input_data:
 for line in input_data:
 interactions_list.append(line.strip('\n').split(', '))
#Get the types of energy reported in the first line
#assuming from PIEDA_mat and all report same energy types
raw_energy_types = interactions_list[0][3]
energy_types = raw_energy_types.split('+')
print energy_types
Etot_index, Ees_index, Eex_index, Edisp_index, Ect_index, Gsol_index = (None,)*6
#Get whatever the ordering is in this list
for i in range(len(energy_types)):
 if energy_types[i] == "Etot":
 Etot_index = i+4
 if energy_types[i] == "Ees":
 Ees_index = i+4
 if energy_types[i] == "Eex":
 Eex_index = i+4
 if energy_types[i] == "Edisp":
 Edisp_index = i+4
 if energy_types[i] == "Ect":
 Ect_index = i+4
 if energy_types[i] == "Gsol":
 Gsol_index = i+4
complexes = []
if Etot_index is not None:
 Etot_engs_list = []
if Ees_index is not None:
 Ees_engs_list = []
if Eex_index is not None:
 Eex_engs_list = []
if Edisp_index is not None:
 Edisp_engs_list = []
if Ect_index is not None:
 Ect_engs_list = []
if Gsol_index is not None:
 Gsol_engs_list = []
#Negative values so bars don't cancel
if Etot_index is not None:
 Etot_engs_neg_list = []
if Ees_index is not None:
 Ees_engs_neg_list = []
if Eex_index is not None:
 Eex_engs_neg_list = []
if Edisp_index is not None:
 Edisp_engs_neg_list = []
 if Ect_index is not None:
 Ect_engs_neg_list = []
if Gsol_index is not None:
 Gsol_engs_neg_list = []
for interaction in interactions_list:

#Terse printing which doesn't give fragment numbers and types
 #complexes.append( interaction[0].split('_')[0] + ' ' + interaction[0].split('_')[1] + ' '
+ interaction[0].split('_')[-1])
#Setting needed for 1bps 1bps differnce
 #complexes.append( interaction[0].split('_')[0] + ' ' + interaction[0].split('_')[2] + ' '
+ interaction[0].split('_')[-1])
#Verbose printing which does give all the details
 complexes.append( interaction[0].split('_')[0] + ' ' + interaction[0].split('_')[1] + ' '
+ interaction[0].split('_')[2] + ' ' + interaction[1] + ' ' + interaction[2])
 if Etot_index is not None:
 if float(interaction[Etot_index]) > 0:
 Etot_engs_list.append(float(interaction[Etot_index]))
 Etot_engs_neg_list.append(0)
 else:
 Etot_engs_list.append(0)
 Etot_engs_neg_list.append(float(interaction[Etot_index]))
 if Ees_index is not None:
 if float(interaction[Ees_index]) > 0:
 Ees_engs_list.append(float(interaction[Ees_index]))
 Ees_engs_neg_list.append(0)
 else:
 Ees_engs_list.append(0)
 Ees_engs_neg_list.append(float(interaction[Ees_index]))
 if Eex_index is not None:
 if float(interaction[Eex_index]) > 0:
 Eex_engs_list.append(float(interaction[Eex_index]))
 Eex_engs_neg_list.append(0)
 else:
 Eex_engs_list.append(0)
 Eex_engs_neg_list.append(float(interaction[Eex_index]))
 if Edisp_index is not None:
 if float(interaction[Edisp_index]) > 0:
 Edisp_engs_list.append(float(interaction[Edisp_index]))
 Edisp_engs_neg_list.append(0)
 else:
 Edisp_engs_list.append(0)
 Edisp_engs_neg_list.append(float(interaction[Edisp_index]))
 if Ect_index is not None:
 if float(interaction[Ect_index]) > 0:
 Ect_engs_list.append(float(interaction[Ect_index]))
 Ect_engs_neg_list.append(0)
 else:
 Ect_engs_list.append(0)
 Ect_engs_neg_list.append(float(interaction[Ect_index]))
 if Gsol_index is not None:
 if float(interaction[Gsol_index]) > 0:
 Gsol_engs_list.append(float(interaction[Gsol_index]))
 Gsol_engs_neg_list.append(0)
 else:
 Gsol_engs_list.append(0)
 Gsol_engs_neg_list.append(float(interaction[Gsol_index]))
#Convert lists to arrays because that's how matplotlib likes it and allows stacked columns
if Etot_index is not None:
 Etot_engs = np.array( Etot_engs_list )
if Ees_index is not None:
 Ees_engs = np.array( Ees_engs_list )
if Eex_index is not None:
 Eex_engs = np.array( Eex_engs_list )
if Edisp_index is not None:
 Edisp_engs = np.array( Edisp_engs_list )
if Ect_index is not None:
 Ect_engs = np.array( Ect_engs_list )
if Gsol_index is not None:
 Gsol_engs = np.array( Gsol_engs_list )
#And their negative counter-parts

if Etot_index is not None:
 Etot_neg_engs = np.array( Etot_engs_neg_list )
if Ees_index is not None:
 Ees_neg_engs = np.array( Ees_engs_neg_list )
if Eex_index is not None:
 Eex_neg_engs = np.array( Eex_engs_neg_list )
if Edisp_index is not None:
 Edisp_neg_engs = np.array( Edisp_engs_neg_list )
if Ect_index is not None:
 Ect_neg_engs = np.array( Ect_engs_neg_list )
if Gsol_index is not None:
 Gsol_neg_engs = np.array( Gsol_engs_neg_list )
x_pos = np.arange(len(complexes))
#Aesthetic parameters
bar_width = 0.75
opacity = 1.0
pos_offset = 0
neg_offset = 0
#Creating the stacked bars
if Etot_index is not None:
 pEtot = plt.bar(x_pos, Etot_engs, bar_width, color='b', alpha=opacity, label='Total
interactions')
 pos_offset += Etot_engs
if Ees_index is not None:
 pEes = plt.bar(x_pos, Ees_engs, bar_width, color='g', alpha=opacity,
label='Electrostatic', bottom=pos_offset)
 pos_offset += Ees_engs
if Eex_index is not None:
 pEex = plt.bar(x_pos, Eex_engs, bar_width, color='r', alpha=opacity, label='Exchangerepulsion', bottom=pos_offset)
 pos_offset += Eex_engs
if Edisp_index is not None:
 pEdisp = plt.bar(x_pos, Edisp_engs, bar_width, color='c', alpha=opacity,
label='Dispersion', bottom=pos_offset )
 pos_offset += Edisp_engs
if Ect_index is not None:
 pEct = plt.bar(x_pos, Ect_engs, bar_width, color='m', alpha=opacity, label='Chargetransfer', bottom=pos_offset )
 pos_offset += Ect_engs
if Gsol_index is not None:
 pGsol = plt.bar(x_pos, Gsol_engs, bar_width, color='y', alpha=opacity, label='Solvation',
bottom=pos_offset )
#And their negative counterparts
if Etot_index is not None:
 pEtot = plt.bar(x_pos, Etot_neg_engs, bar_width, color='b', alpha=opacity)
 neg_offset += Etot_neg_engs
if Ees_index is not None:
 pEes = plt.bar(x_pos, Ees_neg_engs, bar_width, color='g', alpha=opacity,
bottom=neg_offset)
 neg_offset += Ees_neg_engs
if Eex_index is not None:
 pEex = plt.bar(x_pos, Eex_neg_engs, bar_width, color='r', alpha=opacity,
bottom=neg_offset)
 neg_offset += Eex_neg_engs
if Edisp_index is not None:
 pEdisp = plt.bar(x_pos, Edisp_neg_engs, bar_width, color='c', alpha=opacity,
bottom=neg_offset )
 neg_offset += Edisp_neg_engs
if Ect_index is not None:
 pEct = plt.bar(x_pos, Ect_neg_engs, bar_width, color='m', alpha=opacity, bottom=neg_offset
)
 neg_offset += Ect_neg_engs
if Gsol_index is not None:
 pGsol = plt.bar(x_pos, Gsol_neg_engs, bar_width, color='y', alpha=opacity,
bottom=neg_offset )
plt.xticks(x_pos+(bar_width/2), complexes, fontsize=8, rotation='vertical')
plt.ylabel('Energy (kcal/mol)')
plt.legend(loc=8,prop={'size':8})

plt.grid(color='k', linestyle='-', linewidth=0.10)
plt.show()
 

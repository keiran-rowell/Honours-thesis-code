# Written by Keiran Rowell for his Honours project in Chemistry at UNSW.
# Largely written during the second half of 2014.
################################################################################
#This version I added a few checks so free dna sequences are handled when searching
import sys
import os
import re
#gamess_files_dir='/media/21a34721-cb72-4dba-92d3-
27847c51e266/Ashley/facio_files/minimised_bestmember_pdbs/FMO_input/'
gamess_files_dir='/mnt/crashley/Ashley/FMO_calculations/RI-SCS_MP2_6-
31Gd_PCM_PIEDA/output_files/'
gamess_out_files= [f for f in os.listdir(gamess_files_dir) if f.endswith('.out')]
complexes_list = []
#------------------------------------Functions--------------------------------#-
class Fragment(object):
 "fragments from the FMO calculations"
 def __init__(self, frag_type=None, frag_num=None, frag_charge=None, frag_mer=None,
frag_complex=None, frag_job=None, stacked_next=None, stacked_prev=None, frag_paired=None,
stacked_next_leading=None, stacked_next_lagging=None, stacked_prev_leading=None,
stacked_prev_lagging=None):
 self.frag_name = str(frag_num) + '-' + frag_type
 self.frag_type = frag_type
 self.frag_num = frag_num
 self.frag_charge = frag_charge
 self.frag_mer = frag_mer
 self.frag_complex = frag_complex
 self.frag_job = frag_job
 self.stacked_next = stacked_next
 self.stacked_prev = stacked_prev
 self.frag_paired = frag_paired
 #Leading and lagging used for the Chromos since they stack with 4 bases.
 self.stacked_next_leading = stacked_next_leading
 self.stacked_next_lagging = stacked_next_lagging
 self.stacked_prev_leading = stacked_prev_leading
 self.stacked_prev_lagging = stacked_prev_lagging
 
def make_frag(frag_type, frag_num, frag_charge, frag_mer, frag_complex, frag_job):
 fragment = Fragment(frag_type, frag_num, frag_charge, frag_mer, frag_complex, frag_job)
 return fragment
 
def get_complex_details(file):
 complex = {}
 job_name = file
 complex['job_name'] = job_name
 split_job_name = job_name.split("_")
 #assign details to complex's dict
 if split_job_name[1] == 'dual':
 complex['inter_seq'] = split_job_name[0]
 complex['ligand'] = split_job_name[2]
 complex['groove'] = split_job_name[3]
 complex['lig_amount'] = 2 #not sure if this is the right way to implement this check
 elif split_job_name[1] == 'bestmin':
 complex['inter_seq'] = split_job_name[0]
 complex['ligand'] = 'free_dna'
 complex['groove']= None
 complex['lig_amount']= 0
else:
 complex['inter_seq']= split_job_name[0]
 complex['ligand']= split_job_name[1]
 complex['groove']= split_job_name[2]
 complex['lig_amount']= 1
 if complex['ligand'] == 'free_dna':
 complex['complex_name']= complex['inter_seq']+ '_'+ complex['ligand']
else:
 complex['complex_name']= complex['inter_seq']+ '_'+ complex['ligand']+ '_' + complex['groove']
 #verify that it's a valid DNA sequence
if not re.match("^[ACTG]+$", complex['inter_seq']):
 print("Bad complex name, make sure the first field contains only valid DNA characters(A,T,C,G)")
 quit()
if len(complex['inter_seq']) == 4:
 complex['mer']= 14
elif len(complex['inter_seq']) == 3:
 complex['mer']= 13
else:
 print
("Sequence not 13-mer or 14-mer, I don't know how to handle this")
 quit() 
 
 return complex
 
def build_lagging_strand(complex):
lagging_strand_seq= ""
leading_strand_seq= complex['leading']

for base in leading_strand_seq:
 if base == "A":
  lagging_strand_seq += "T"
 elif base == "C":
  lagging_strand_seq += "G"
 elif base == "T":
  lagging_strand_seq += "A"
 elif base == "G":
  lagging_strand_seq += "C"
 
 complex['lagging'] = lagging_strand_seq

 return complex

def build_duplex(complex):
 reversed_lagging = complex['lagging'][::-1]
 complex['duplex_seq'] = complex['leading'] + reversed_lagging
 
 return complex
 
def determine_number_sodiums(complex):
 if complex['ligand'] == 'C3NC3': #All other arrangements lead to two charges, this one three.
  complex['num_sodiums']= complex['mer'] * 2 - 5
 elif complex['ligand'] == 'free_dna':
  complex['num_sodiums']= complex['mer'] * 2- 2
 else:
  complex['num_sodiums']= complex['mer']* 2 - 4
  
 return complex
 
def determine_number_ligand_frags(complex):
 if complex['lig_amount'] == 2:
  complex['num_lig_frags'] = 2
 elif complex['ligand'] == 'free_dna':
  complex['num_lig_frags'] = 0
 else:
  complex['num_lig_frags'] = 3 #single ligands in my case are bisintercalators, broken up into three fragments.
  
 return complex
def get_frag_charges_and_labels(complex):
 """Takes in the name of a file (job_name) and searches the GAMESS .out for the Fragment
statistics section. Populates a list with the fragment labels (complex['frag_labels']) and
charges (complex['frag_charges']) present in the .out file, rather than relying on searching
for the matching .inp file, or it being echoed."""
 job_name = complex['job_name']

 raw_frag_stats_list = []
 clean_frag_stats_list = []
 target_file = gamess_files_dir + job_name

 with open(target_file) as input_data:
 for line in input_data:
 if line.startswith(' I NAME Q NAT0 NATB NA NAO LAY MUL SCFTYP
NOP MOL CONV'): #Line just after Fragment statistics header
 break
 for line in input_data:
 if line.startswith(' Close fragment pairs, distance relative to vdW radii'):
 break
 raw_frag_stats_list.append(line.replace('=','').split()) #Remove = seperator
 input_data.close()
 clean_frag_stats_list = filter(None, raw_frag_stats_list)

 frag_labels_list = []
 frag_charges_list = []
 for i in range(len(clean_frag_stats_list)):
 frag_labels_list.append(clean_frag_stats_list[i][1])
 frag_charges_list.append(clean_frag_stats_list[i][2])
 complex['frag_labels'] = frag_labels_list
 complex['frag_charges'] = frag_charges_list
 return complex
def populate_fragment_list(complex):
 complex['fragments'] = []
 #Need to grab chromophore locations from labels list
 #My definition is CR1 is the first chromo encounter moving along the 'tip' (CGATG)
 for j in range(len(complex['frag_labels'])):
 if complex['frag_labels'][j] == 'CR1':
 CR1_index = j
 if complex['frag_labels'][j] == 'CR2':
 CR2_index = j
 #If I have forgotten to label CR1 and CR2 this will cause errors, but that will point out
my mistake. 

#Populate fragment list of complex with fragment objects, should make this a separate
function
 for i in range( (2*len(complex['duplex_seq'])) + complex['num_lig_frags'] +
complex['num_sodiums']):
 #print complex['frag_charges'][i]
 if i < (2*len(complex['duplex_seq'])):
 if i % 2 == 0:
 new_fragment = make_frag('p', i+1, complex['frag_charges'][i], complex['mer'],
complex['complex_name'], complex['job_name'])
 else:
 if complex['duplex_seq'][i/2] == "A":
 new_fragment = make_frag('A', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 elif complex['duplex_seq'][i/2] == "C":
 new_fragment = make_frag('C', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 elif complex['duplex_seq'][i/2] == "T":
 new_fragment = make_frag('T', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 elif complex['duplex_seq'][i/2] == "G":
 new_fragment = make_frag('G', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 elif i >= (2*len(complex['duplex_seq'])) and i < ( (2*len(complex['duplex_seq'])) +
complex['num_lig_frags'] ):
 if i == CR1_index:
 new_fragment = make_frag('CR1', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 elif i == CR2_index:
new_fragment = make_frag('CR2', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 else:
 new_fragment = make_frag('LNK', i+1, complex['frag_charges'][i],
complex['mer'], complex['complex_name'], complex['job_name'])
 else:
 new_fragment = make_frag('Na', i+1, complex['frag_charges'][i], complex['mer'],
complex['complex_name'], complex['job_name'])
 complex['fragments'].append(new_fragment)

 #Sanity check to make sure you've created the right number of units.
 if len(complex['fragments']) != len(complex['frag_charges']):
 print("Something went very wrong, the number of fragment objects doesn't match the
.out file's charge list!")
 quit()
 return complex
def get_stacked(complex):
 #Step through to find the chromophore fragments
 for fragment in complex['fragments']:
 if fragment.frag_type == 'CR1':
 CR1_frag = fragment
 if fragment.frag_type == 'CR2':
 CR2_frag = fragment
 
if complex['ligand'] != 'free_dna':
 if complex['mer'] == 14 and complex['ligand']:
 for i in range(1,56,2):
 #First define all the end fragment sections
 if i == 1:
 complex['fragments'][i].stacked_prev = None
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 27:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = None
 elif i == 29:
 complex['fragments'][i].stacked_prev = None
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 55:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = None
 #Then rules for the bases stacked with the chromophores
 elif i == 11:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR1_frag
 elif i == 13:
 complex['fragments'][i].stacked_prev = CR1_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 15:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR2_frag
 elif i == 17:
 complex['fragments'][i].stacked_prev = CR2_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 39:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR2_frag
 elif i == 41:
 complex['fragments'][i].stacked_prev = CR2_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 43:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR1_frag
 elif i == 45:
 complex['fragments'][i].stacked_prev = CR1_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 #Finally the 'standard' inter-strand stacks
 else:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 #Deal with the 4 stacking of chromophores by defining major and minor strand
 CR1_frag.stacked_next_leading = complex['fragments'][13]
 CR1_frag.stacked_next_lagging = complex['fragments'][45]
 CR1_frag.stacked_prev_leading = complex['fragments'][11]
 
 CR1_frag.stacked_prev_lagging = complex['fragments'][43]
 CR2_frag.stacked_next_leading = complex['fragments'][17]
 CR2_frag.stacked_next_lagging = complex['fragments'][41]
 CR2_frag.stacked_prev_leading = complex['fragments'][15]
 CR2_frag.stacked_prev_lagging = complex['fragments'][39]

 elif complex['mer'] == 13:
 for i in range(1,52,2):
 #First define all the end fragment sections
 if i == 1:
 complex['fragments'][i].stacked_prev = None
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 25:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = None
 elif i == 27:
 complex['fragments'][i].stacked_prev = None
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 51:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = None
 #Then rules for the bases stacked with the chromophores
elif i == 11:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR1_frag
 elif i == 13:
 complex['fragments'][i].stacked_prev = CR1_frag
 complex['fragments'][i].stacked_next = CR2_frag
 elif i == 15:
 complex['fragments'][i].stacked_prev = CR2_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 elif i == 37:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = CR2_frag
 elif i == 39:
 complex['fragments'][i].stacked_prev = CR2_frag
 complex['fragments'][i].stacked_next = CR1_frag
 elif i == 41:
 complex['fragments'][i].stacked_prev = CR1_frag
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 #Finally the 'standard' inter-strand stacks
 else:
 complex['fragments'][i].stacked_prev = complex['fragments'][i-2]
 complex['fragments'][i].stacked_next = complex['fragments'][i+2]
 #Deal with the 4 stacking of chromophores by defining major and minor strand
 CR1_frag.stacked_next_leading = complex['fragments'][13]
 CR1_frag.stacked_next_lagging = complex['fragments'][41]
 CR1_frag.stacked_prev_leading = complex['fragments'][11]
 CR1_frag.stacked_prev_lagging = complex['fragments'][39]
 CR2_frag.stacked_next_leading = complex['fragments'][15]
 CR2_frag.stacked_next_lagging = complex['fragments'][39]
 CR2_frag.stacked_prev_leading = complex['fragments'][13]
 CR2_frag.stacked_prev_lagging = complex['fragments'][37]
 else:
 print "%s is not a 14- or 13-mer, can't deal" % complex['complex_name']

 return complex
def get_paired(complex):
 #Get the fragment which would correspond to the Watson-Crick paired base
 if complex['mer'] == 14:
 for i in range(1,56,2):
 complex['fragments'][i].frag_paired = complex['fragments'][56-i]
 elif complex['mer'] == 13:
 for i in range(1,52,2):
 complex['fragments'][i].frag_paired = complex['fragments'][52-i]
 else:
 print "%s is not a 14- or 13-mer, can't deal" % complex['complex_name']
 return complex
def get_twobody_matrix(complex):
 """Takes in the name of a file (job_name) and searches the GAMESS .out for the final 'Twobody FMO properties' section. Returns a 2D array (list) of the interaction energies."""
 job_name = complex['job_name'] 
 
  twobody_matrix = []
 twobody_header_count = 0
 target_file = gamess_files_dir + job_name
 with open(target_file) as input_data:
 for line in input_data:
 if line.startswith(' I J DL Z R Q(I->J) EIJ-EI-EJ dDIJ*VIJ total
Ees Eex Ect+mix Edisp Gsol'):
 twobody_header_count += 1
 if twobody_header_count > 1 : #Hacky currenly, would like to do reverse
search from bottom of file but difficult with incrementers
 input_data.next() #Need to remove the ---- header, just step silently
along one more iteration
 break
 for line in input_data:
 if line.startswith('\n'): #Continue reading until the newline space at the end
of the matrix
 break
 twobody_matrix.append(line.split()) #Concise way to add all these lines to the
matrix
 complex['matrix'] = twobody_matrix

 return complex
def get_IFIE(frag_a, frag_b, energy_request, complexes_list):
 """ Given two fragments follwed by an arguement of one or more energy types (Ees+Eex,
Etot, etc.),
 will return a list of the interaction energy between the two in the format:
 [COMPLEX_NAME, FRAG_A, FRAG_B, ENERGYTYPES, ENGVAL1, ENGVAL2, ..., ENGVALN] """
 #Sanity check that calling two fragments in the same complex
 if frag_a.frag_job != frag_b.frag_job:
 print("Can't get IFIE on these fragments %s, %s: not from the same job!" ) %
(frag_a.frag_name, frag_b.frag_name)#Should probably be an exception
 return None

 #Split up requested energies so that the matrix can be searched.
 queried_energies = energy_request.split('+')
 energy_positions_dict = {'Etot':8, 'Ees':9, 'Eex':10, 'Ect':11, 'Edisp':12, 'Gsol':13}
 energy_pos_vals = []
 for energy in queried_energies:
 energy_pos_vals.append(energy_positions_dict[energy])
 IFIE_line = []
 for complex in complexes_list:
 if complex['job_name'] == frag_a.frag_job:
 target_complex = complex
 target_complex_matrix = target_complex['matrix']
 for line in target_complex_matrix:
 if (eval(line[0]) == frag_a.frag_num and eval(line[1]) == frag_b.frag_num) or
(eval(line[0]) == frag_b.frag_num and eval(line[1]) == frag_a.frag_num): #I think I need to
clean this up by looking at the ordering
 IFIE_line = [frag_a.frag_complex, frag_a.frag_name, frag_b.frag_name,
energy_request]
 for position in energy_pos_vals:
 IFIE_line.append(line[position])
 return IFIE_line
#Should add a toggle hear to give full frag details or not
def print_IFIE(frag_a, frag_b, energy_request, complexes_list):
 IFIE_line = get_IFIE(frag_a, frag_b, energy_request, complexes_list)
 print str(IFIE_line).translate(None,"'").strip("[").strip("]") #remove the list
demarcations for easier reading/parsing
#---------------------------------------------
Main---------------------------------------------#
#I have even more separating out into functions to do.
for file in gamess_out_files:
 complex = get_complex_details(file)
# print("Adding %s ...") % complex['job_name']
 complex = build_leading_strand(complex)
 complex = build_lagging_strand(complex)
 complex = build_duplex(complex)
 complex = determine_number_sodiums(complex)
 complex = determine_number_ligand_frags(complex)
 complex = get_frag_charges_and_labels(complex)
# print complex #before you get the horrendous amount of numbers from the matrix
 complex = populate_fragment_list(complex)
 complex = get_stacked(complex)
 
  complex = get_paired(complex)
 complex = get_twobody_matrix(complex)
 complexes_list.append(complex)

#print('\n')
#print("The complexes loaded in are:")
#for complex in complexes_list:
# print complex['job_name']
#-----------------------------------------Searching
logic----------------------------------------#
#This way I can just write seperate pieces of logic and just supply them as an argument
logic_code = sys.stdin.read()
exec logic_code
 

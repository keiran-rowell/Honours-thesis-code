#!/bin/bash
#WIP. This script will take a frame from clustering analysis and submitting to a brief MM
minimisation (steepest descent) in order remove transient unfavourable interactions.
#Alongside minimising the bestmember it also handles: stripping solvation, adding counter-ions
and converting naming scheme to ready it for facio
cd ./clustering
pop_clust=$(sort -n -r -k3 Centroids_3A.stats | head -n 1 | awk '{print $1}')
echo "The most populous cluster is $pop_clust"
abs_frame=$(sort -n -k2 centroid0${pop_clust}.member.dat | head -n 1 | awk '{print $1}')
echo "The besmember frame is at frame number $abs_frame"
#200 frames be mdcrd and the first 22 are equilibriation
div_num=$(($abs_frame/200))
crd_file_num=$(($div_num + 23))
echo "It is in the original mdcrd file number $crd_file_num"
mod_num=$(($abs_frame%200))
frame_num=$mod_num
echo "In that mdcrd file it is frame $mod_num"
echo "Extracting frame..."
cd ..
SCRIPTS_DIR=/home/keiran/scripts
CONTAINING_DIR=${PWD##*/}
NAME=$CONTAINING_DIR
sed -e s/CRDNUM/"$crd_file_num"/g -e s/FRAMENUM/"$frame_num"/g -e s/NAME/"$NAME"/g < $
{SCRIPTS_DIR}/extract_frame_from_clust_template.trajin > $
{SCRIPTS_DIR}/extract_frame_from_clust_edited.trajin
#Check whether gziped or not. Under construction
if [[ ! -f md${crd_file_num}.mdcrd.gz ]]
then
 sed -i s/.gz//g ${SCRIPTS_DIR}/extract_frame_from_clust_edited.trajin
 if [[ ! -f md${crd_file_num}.mdcrd ]]
 then
 echo "The coordinate file md${crd_file_num}.mdcrd appears to be missing"
 exit
 fi
fi
/usr/local/amber12/bin/ptraj ${NAME}_wat.prmtop < $
{SCRIPTS_DIR}/extract_frame_from_clust_edited.trajin
mkdir ./clustering/minimisation
cp "extracted_${NAME}_frame_${crd_file_num}_${frame_num}.restrt.${frame_num}"
./clustering/minimisation
#have to put the frame number after the extension due to quirk in ptraj
cp ${NAME}_wat.prmtop ./clustering/minimisation
cd ./clustering/minimisation
cp "extracted_${NAME}_frame_${crd_file_num}_${frame_num}.restrt.${frame_num}" "min0.restrt"
cp ${SCRIPTS_DIR}/min_steepdesc_50cyc.in .
cp ${SCRIPTS_DIR}/repeat_min .
sed -i s/NAME/"${NAME}_wat"/g repeat_min
num_repeats=5
./repeat_min $(($num_repeats -1 )) #Andre's repeat scripts all use a gt comparison
ambpdb -p ${NAME}_wat.prmtop < min${num_repeats}.restrt > ${NAME}_bestmember_minimised.pdb
cp ${SCRIPTS_DIR}/trajin_strip_solvation.in .
sed -i s/SOLVATEDFILE/"${NAME}_bestmember_minimised"/g trajin_strip_solvation.in
sed -i s/STRIPPEDFILE/"${NAME}_bestmember_minimised_stripped"/g trajin_strip_solvation.in
cpptraj ${NAME}_wat.prmtop < trajin_strip_solvation.in
FMO_facio_input_conversion.sh ${NAME}_bestmember_minimised_stripped.mol2
crashley=`echo $crashley`
PDB_DIR=${crashley}/facio_files/minimised_bestmember_pdbs
cp ${NAME}_bestmember_minimised_stripped_facio_ready.pdb $PDB_DIR
#should append a label onto my minimised pdbs to say how many cycles they have been minimised for. 

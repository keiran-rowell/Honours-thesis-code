#!/bin/bash
SCRIPTS_DIR=/home/keiran/scripts
units_dir=/home/keiran/units
file="$1"
basename=${file%.mol2}
#Variables for name sanity checking
type=unknown
mer=unkown
sequence=unknown
unit=unkown
groove=unknown
num_fields=$(echo "$file" | awk -F '_' '{print NF}' )
if [ $num_fields -eq 4 ]
 then
 echo "Likely is free DNA"
 type=free_dna
 sequence=$(echo "$file" | awk -F '_' '{print $1}' )
 unit=none
 groove=none
elif [ $num_fields -eq 6 ]
 then
 echo "DNA with drug"
 type=dna_drug
 sequence=$(echo "$file" | awk -F '_' '{print $1}' )
 unit=$(echo "$file" | awk -F '_' '{print $2}' )
 groove=$(echo "$file" | awk -F '_' '{print $3}' )
elif [ $num_fields -eq 7 ]
 then
 if [[ ! "$(echo "$file" | awk -F '_' '{print $1}' )" -eq "dual" ]]
 then
 type=invalid
 echo "If the complex is dual intercalated, 2nd field separated by underscore should
be "dual". Otherwise file naming issue"
 exit
 fi
 echo "DNA with 2 monointercalators"
 type=dna_2mono
 sequence=$(echo "$file" | awk -F '_' '{print $1}' )
  unit=$(echo "$file" | awk -F '_' '{print $3}' )
 groove=$(echo "$file" | awk -F '_' '{print $4}' )
else
 echo "I'm confused! The input file isn't named correctly, wrong amount of fields separated
by underscores."
 type=invalid
 exit
fi
#Regex to make sure valid dna sequence
if [[ $sequence =~ ^[ACTG]+$ ]]
then
 echo "Valid sequence"
else
 type=invalid
 echo "Invalid sequence, please make sure file name is correct"
 exit
fi
if [[ ${#sequence} -eq 4 ]]
 then mer=14
elif [[ ${#sequence} -eq 3 ]]
 then mer=13
else
 type=invalid
fi
echo "type is $type"
echo "sequence name is $sequence"
echo "unit name is $unit"
echo "groove side is $groove"
echo "mer is ${mer}-mer"
#Arithmetic to determine the number of counter-ions to add
if [[ $mer -eq 14 ]]
 then
 num_ions=26
elif [[ $mer -eq 13 ]]
 then
 num_ions=24
else
 type=invalid
 echo "This seems to be the wrong type of -mer! Please check file name."
 exit
fi
if [ "$type" == "free_dna" ];
 then
 :
elif [ "$type" == "dna_drug" ];
 then
 if [ "$unit" == "C3NC3" ] #Check that it's not the +3 ligand rather than all the other +2
ligands
 then
 num_ions=$(($num_ions-3))
 else
 num_ions=$(($num_ions-2))
 fi
elif [ "$type" == "dna_2mono" ];
 then
 num_ions=$(($num_ions-2))
else
 type=invalid
 echo "Problem getting complex type, please check naming conventions."
fi
if [ "$type" == "invalid" ];
 then
 echo "Somewhere along the way the complex type became invalid, check your naming
conventions."
fi
if [ "$unit" == "C3NC3" ]
 then
 unit="C3" #Sadly the Amber programs seem to cap at 3 letter unit
fi
if [ "$unit" == "9aa" ]
 then
 unit="9a"
fi
if [ "$unit" == "9AA" ]
 then
 unit="9a"
fi
#See what the program has decided
echo "Complex -mer is ${mer}-mer"
echo "Complex type is $type"
echo "Therefore adding $num_ions Na+ counter-ions"
#Fix hydrogen naming incopmatibilites and replace MOL with unit type
cp $file ${basename}_fixed.mol2
#Below was need when I was using .pdb files instead of .mol2 files
#sed -e s/MOL/"${unit}L"/g -e s/\'HO5/HO5\'/g -e s/\'HO3/HO3\'/g -e s/\'H5\'/H5\'\'/g -e
s/\'H2\'/H2\'\'/g < "$file" > "${basename}_fixed.mol2"
#Create leap template for adding ions in tleap, then add those counter-ions
sed -e s#UNITS_DIR#"$units_dir"#g -e s/UNIT/"${unit}L"/g -e s/BASENAME/"$basename"/g -e
s/NUMIONS/$num_ions/g < ${SCRIPTS_DIR}/leap_addions_template.cmd > $
{SCRIPTS_DIR}/leap_addions_edited.cmd
tleap -f ${SCRIPTS_DIR}/leap_addions_edited.cmd
#Fix more naming incompatibles with facio
sed -e s/HO5\'/\ H5\'/g -e s/HO3\'/\ H3\'/g -e s/\'H5/\ H5/g -e s/\'H2/\ H2/g -e s/\ Pt/Pt\ /g
< "${basename}_ions.pdb" > "${basename}_facio_ready.pdb"

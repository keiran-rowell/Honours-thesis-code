#!/bin/bash
#Automating the clustering procedure as outline by Ross Walker
#(http://ambermd.org/tutorials/basic/tutorial3/section6.htm)
#Setting up variables:
#3 works best empirically for my system.
RADIUS_VALUE=3
#truncate full path to current directory
CONTAINING_DIR=${PWD##*/}
#usually the .prmtop and dir follow same naming convertion. Alter if not,
NAME=$CONTAINING_DIR
echo "Clustering of trajectory of $NAME"
#Alter to relevant directory with all cluster scripts. Make sure MMTSB is accesible
SCRIPTS_DIR=/home/keiran/scripts
echo "Using scripts from $SCRIPTS_DIR"
#Perform job:
#gunzip required .prmtop files:
if [ -f "${NAME}_wat.prmtop.gz" ]
then
 gunzip ${NAME}_wat.prmtop.gz
fi
if [ -f "${NAME}_vac.prmtop.gz" ]
then
 gunzip ${NAME}_vac.prmtop.gz
fi
#Check for .prmtop
if [ ! -f "${NAME}_wat.prmtop" ]
then
 echo "${NAME}_wat.prmtop does not exist or is misnamed"
 exit
fi
if [ ! -f "${NAME}_vac.prmtop" ]
then
 echo "${NAME}_vac.prmtop does not exist or is misnamed"
 exit
fi
#Create the binops file to work with later.
#Should add exit status if no md23.mdrcrd(.gz)
if [ -f "md23.mdcrd.gz" ]
then
 cpptraj ${NAME}_wat.prmtop < ${SCRIPTS_DIR}/mdcrdgz_to_binpos.ptraj > mdcrd_to_binpos.out
elif [ -f "md23.mdcrd" ]
then
 cpptraj ${NAME}_wat.prmtop < ${SCRIPTS_DIR}/mdcrd_to_binpos.ptraj > mdcrd_to_binpos.out
else
 echo "No md23.mdcrd file! Aborting."
 exit
fi
#make directories to work in.
mkdir clustering
cd clustering
mkdir PDBfit
#generate PDBs for each frame (20,000 for 10ns)
/usr/local/amber12/bin/ptraj ../${NAME}_vac.prmtop < ${SCRIPTS_DIR}/extract_pdb.ptraj
#sane numbering, adding leading 0s
cd ./PDBfit
${SCRIPTS_DIR}/fix_numbering_pdb.csh
rm complex_clust.pdb #An un-numbered pdb causes segfault
#Standard k-means clustering
rm ../clustfils
ls -1 . > ../clustfils
kclust -mode rmsd -centroid -cdist -heavy -lsqfit -radius ${RADIUS_VALUE} -maxerr 1
-iterate ../clustfils > ../Centroids_3A
#extract centroids
cd ..
awk -f ${SCRIPTS_DIR}/extract_centroids.awk Centroids_3A | tee Centroids_3A.stats

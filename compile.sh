#!/bin/bash

echo compiling in $STAR_VERSION
rsync -va $STAR/StRoot/StEvent ./StRoot/
rsync -va $STAR/StRoot/StMuDSTMaker ./StRoot/

#remove the kfp flag - old way
#sed -i 's/#define\s\+__kfpAtFirstHit__//g' StRoot/StEvent/StTrackTopologyMap.h

#instead patch /afs/rhic.bnl.gov/star/packages/SL21a from official production - also solving topoMap problem
cp /afs/rhic.bnl.gov/star/packages/SL21a/StRoot/StEvent/StTrackTopologyMap.* ./StRoot/StEvent/
sed -i 's/#define kiTpcIdentifier.*/#define kiTpcIdentifier 39/g' StRoot/StEvent/StDetectorDefinitions.h  

cons


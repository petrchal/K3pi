#!/bin/bash

echo chp0
codeDir=`pwd`
cd $HOME/kaons_resKFP2
productionId=`date +%F_%H-%M`_$1
analyzer="petrchal"

mkdir $productionId
cd $productionId

#copy needed folders
cp -r $codeDir/.sl73* ./
cp -Lr $codeDir/StRoot ./
cp $codeDir/StRoot/macros/kfpAnalysis.C ./
cp $codeDir/StRoot/macros/lMuDst.C ./
cp $codeDir/setDEV2.csh ./
mkdir starSubmit
cp $codeDir/SimAnalysis.xml ./starSubmit

mkdir -p production
mkdir -p report
mkdir -p csh
mkdir -p list
mkdir -p jobs
mkdir -p jobs/log
mkdir -p jobs/err

path=`pwd -P`
path=$( echo $path | sed 's|//|/|g' )

baseFolder=${path}
echo $baseFolder

noPID=true

star-submit-template -template ./starSubmit/SimAnalysis.xml -entities basePath=${baseFolder},prodId=${productionId},noPID=${noPID}

#!/bin/bash
# Dave Hill, OPenPedCan 2023
# Run OpenPedCan modules to generate MTP tables

set -e
set -o pipefail

printf "Start generating mtp tables...\n\n"

#Set the desired folder locations for each file type
csv_diseases_folder="/home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-csv/diseases"
csv_targets_folder="/home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-csv/targets"
json_diseases_folder="/home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-json/diseases"
json_targets_folder="/home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-json/targets"


#Recursively download all of the diseases json files into the mtp-json/diseases scratch folder
wget -P /home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-json/ --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/diseases


#Recursively download all of the targets json files into the mtp-json/targets scratch folder
wget -P /home/ubuntu/volume/OpenPedCan-analysis/scratch/mtp-json/ --recursive --no-parent --no-host-directories --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/22.11/output/etl/json/targets

#Move to the mtp-json/diseases folder
cd $json_diseases_folder

#Make an array of the names of each json file in the mtp-json/diseases folder
dirlist=(`ls ${prefix}*.json`)


#Runs a loop to change each mtp-json/diseases/ file into a csv file
for filelist in $json_diseases_folder/*.json
do 
	echo "$(basename  $filelist .json)"
	name=$(basename  $filelist .json)
	echo "$filelist"
	dasel -r json -w csv < $filelist > $name.csv

done

#Move all csv files in the mtp-json/diseaes folder to the mtp-csv/diseases folder
mv $json_diseases_folder/*.csv $csv_diseases_folder/.

#Move to the mtp-json/targets folder
cd $json_targets_folder

##Make an array of the names of each json file in the mtp-json/targets folder
dirlist=(`ls ${prefix}*.json`)

#Runs a loop to change each mtp-json/targets/ file into a csv file
for filelist in $json_targets_folder/*.json
do
        echo "$(basename  $filelist .json)"
        name=$(basename  $filelist .json)
        echo "$filelist"
        dasel -r json -w csv < $filelist > $name.csv

done

#Move all csv files in the mtp-json/targets folder to the mtp-csv/targets folder
mv $json_targets_folder/*.csv $csv_targets_folder/.

Rscript ~/volume/OpenPedCan-analysis/analyses/mtp-annotations/01-mtp-annotations.R 

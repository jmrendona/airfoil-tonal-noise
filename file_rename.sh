#!/bin/bash

master_directory="/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/Others/Mine/2025-ISAE-Exp/NACA0015/Cp/Exp-NACA0015-2/"

cd $master_directory

for i in {1..47};do #Iteration in mics
	mv "NACA0015-SVA_sweep-"$i".csv" "NACA0015_SVA_sweep-"$i".csv"

done

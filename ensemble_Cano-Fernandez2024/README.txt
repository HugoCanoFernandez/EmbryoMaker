

This is the code to reproduce the ensemble of 23,432 initial conditions presented in the article Cano-Fernandez et al. 2024

Requirements:
- Linux operative system (tested for Ubuntu 24.04)
- the compiler gfortan must be installed

To reconstruct the ensemble:

1. Extract the directories ics.tar.xz and new_ics.tar.xz

2. give permissions to the script create_ensemble.sh with the command line:
chmod 777 create_ensemble.sh

3. Then just run the script create_ensemble.sh with:
./create_ensemble.sh


All the 23,432 initial conditions will be created. It may take a while...
The files will be idenfified with a code like the following

icoepi_nuclei.datdyn00001_02256_00002_00006.dat
 
where: 00001 refers to the number of the initial condition (1 to 23,432)
       02256 is the identifier of the network of gene-gene interactions
       00002 and 00006 refer to the cell properties and behaviours that are activated. For example 00002 means no cell property or behavior is activated by one of the genes with a complex pattern of gene expression and 00006 means that the other gene is regulating the adhesion radius of cells. The complete list of cell properties and behaviors and their codes can be found in the file ranges_ENSEMBLE_epi.dat.
       
The initial condition files must be run with EmbryoMaker to obtain the final morphologies

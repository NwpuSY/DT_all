This program is for four datasets : E_data.mat, IC_data.mat, GPCR_data.mat and NR_data.mat
Main Function : Predicting_DT_main.m

Each mat file contains the following variables:

DTI --> drug-target interaction matrix
d_s --> chemical-structure-based similarity matrix
d_ATC--> ATC-based similarity matrix
t_s --> sequence-alignment-based similarity matrix
t_pathway_sim-->pathway-based similarity matrix
t_go_sim-->GO-based similarity matrix
t_s_Class--> FC-based similarity matrix
DrugName--> all names of drugs corresponding to the rows of DTI
TargetName--> all names of targets corresponding to the columns of DTI

NOTE: Users should run the main function in MATLAB . 
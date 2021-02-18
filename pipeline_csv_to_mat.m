%Pipe Line convert .csv to .mat 
%by jhon vargas
clc;close all;clear all;
M=readtable(['Base de datos\spotify_pro_5.csv']);
data = table2array(M);
save('spotify_pro_5.mat','data');
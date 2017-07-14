clear all;
close all;
addpath(genpath('lib'))
% model name
name='p53';
% crating model equations based on model equations; should be uncommented after used once 
create(name);
addpath(genpath(['models/',name]))

% initial condition 
y0(1)=0;
y0(2)=0;
y0(3)=0;
y0(4)=1;
y0(5)=1;
y0(6)=1;
y0(7)=0;

y0(8)=0;
y0(9)=0;



% plotting model trajectories. Parametrs are as follows: model name, number
% of observations, time distance between observations, time of initial
% observation,initial conditions, vector with indices of observed variables

plottraj(name,20,1,8,y0,[1, 2, 3])


% Calculating Fisher Information Matrices (FIMs). Parametrs are as follows: model name, number
% of observations, time distance between observations, time of initial
% observation,initial conditions, vector with indices of observed
% variables, measurement error variance, 'All@ for all three data types,
% 'TRUE' for logarithmic parametrisation

F=Fisher(name,20,1,8,y0,[1, 2, 3],10,'All','T');
% extracting the FIMs for three data types 
FIM_TS=F{1};
FIM_TP=F{2};
FIM_DT=F{3};


% Contour plot of the FIM for time-series (TS) and deterministic model (DT) for
% parameters 1 and 2
figure
i=1;
j=2;
plotdFMcc(name, 1,1,FIM_TS,FIM_DT,i,j,'TS','DT'); 

% The same as above but as a heatmap  
figure
plotdFM(name, 1,1,FIM_TS,FIM_DT,i,j);


% FIMs for all parameters as heatmap

plotdFMAll(name, 1,1,FIM_TS,FIM_DT); % First element plotted (TS) polted as heatmap, second (DT) as contourplot



% The same as above but as a contour plot
plotdFMAllcc(name, 1,1,FIM_TS,FIM_DT,'TS','DT')


% sensitivity coeffitients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating sensitifity coeffitients for based on FIM for TS data 
SC1=sensitivities(FIM_TS);

% Plotting sensitifity coeffitients for all data types 
sensitivitiesAll(name, F);

% plotting Cramer-Rao bounds for all parameters for TS and time-point data (TP) 
diaginv(name,F);


% plotting contributions of individual parameters to each of the eigenvalues 
decomp(name,F);











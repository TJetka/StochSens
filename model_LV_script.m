clear all;
close all;
addpath(genpath('lib'))
% model name
name='LV';
% crating model equations based on model equations; should be uncommented after used once 
create(name);
addpath(genpath(['models/',name]))

% initial condition 
y0(1)=1000; % mean of variable 1
y0(2)=1000; % mean of variable 2
y0(3)=1;    % variance of variable 1
y0(4)=1;    % variance of variable 2
y0(5)=0;    % covvariance of variable 1 and 2


N=15; % number of observations
freq=0.1; % time distance between observations
init_T=2; % initial time
obsv=[1, 2]; % indices of observed variables
sigma=10; %  variance of measurement error

% plotting model trajectories. Parametrs are as follows: model name, number
% of observations, time distance between observations, time of initial
% observation,initial conditions, vector with indices of observed variables
plottraj(name,N,freq,init_T,y0,obsv)



% Calculating Fisher Information Matrices (FIMs). Parametrs are as follows: model name, number
% of observations, time distance between observations, time of initial
% observation,initial conditions, vector with indices of observed
% variables, measurement error variance, 'All@ for all three data types,
% 'TRUE' for logarithmic parametrisation

F=Fisher(name,N,freq,init_T,y0,[1, 2],sigma,'All','TRUE');

% extracting the FIMs for three data types 
FIM_TS=F{1};
FIM_TP=F{2};
FIM_DT=F{3};


% Contour plot of the FIM for time-series (TS) and deterministic model (DT) for
% parameters 1 and 2


rx = 1;  % maximal x argument for a quadratic form 
ry = 2;  % maximal y argument for a quadratic form
i = 3; % index of a first parameter
j = 1; % index of a second parameter
FIM_TS = F{1};
FIM_DT = F{2};
figure
plotdFMcc(name, rx, ry, FIM_TS, FIM_DT, i, j, 'TS', 'DT'); 


% The same as above but as a heatmap  
figure
plotdFM(name, rx,ry,FIM_TS,FIM_DT,i,j);


% FIMs for all parameters as heatmap

plotdFMAll(name, rx,ry,FIM_TS,FIM_TP); % First element plotted (TS) polted as heatmap, second (DT) as contourplot



% The same as above but as a contour plot
plotdFMAllcc(name, rx,ry,FIM_TS,FIM_DT,'TS','DT')


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
addpath(genpath('lib'))
name='JAKSTAT';
create(name);

addpath(genpath(['models/',name]))
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters

% converting model paramars into molecular numbers units
W_C=6.02*(10^23)*1400*10^(-18)*10^(-6); % conversion constant (from mM to molecular numbers)
par(2)=par(2)/W_C;

% defining  times to calculate contributions
size_var=14;
y0=zeros(1,(size_var^2-size_var)*0.5 + 2*size_var );
y0(1)=1.99989483593865340000e+002*W_C; % litrature value from Swameye et al. 2003

plottraj(name,50,1,0,y0,[14]);

F=Fisher(name,40,1,10,y0,[14],10,'TS','FALSE');
FIM_TS=F;
F=Fisher(name,40,1,10,y0,[14],10,'TP','FALSE');
FIM_TP=F;
F=Fisher(name,40,1,10,y0,[14],10,'DT','FALSE');
FIM_DT=F;



plotdFMAll(name, 1,1,FIM_TS,FIM_DT); % First element plotted (TS) polted as heatmap, second (DT) as contourplot


i=1;
j=3;
plotdFMcc(name, 1,1,FIM_TS,FIM_DT,i,j,'TS','DT'); 
plotdFM(name, 1,1,FIM_TS,FIM_DT,i,j);



plotdFMAllcc(name, 1,1,FIM_TS,FIM_DT,'TS','DT')


% sensitivity coeffitients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=Fisher(name,40,1,10,y0,[14],10,'All','FALSE');
SC1=sensitivities(F{1});
sensitivitiesAll(name, F);
diaginv(name,F);
decomp(name,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

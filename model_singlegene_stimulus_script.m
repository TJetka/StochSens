clear all;
close all;
addpath(genpath('lib'))
name='singlegene_stimulus';
create(name);

addpath(genpath(['models/',name]))
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters

% defining  times to calculate contributions
size_var=2;
y0=zeros(1,(size_var^2-size_var)*0.5 + 2*size_var );

plottraj(name,40,0.1,0,y0,[2]);

F=Fisher(name,40,0.1,0.1,y0,[2],10,'TS','FALSE');
FIM_TS=F;
F=Fisher(name,40,0.1,0.1,y0,[2],10,'TP','FALSE');
FIM_TP=F;
F=Fisher(name,40,0.1,0.1,y0,[2],10,'DT','FALSE');
FIM_DT=F;



plotdFMAll(name, 1,1,FIM_TS,FIM_DT); % First element plotted (TS) polted as heatmap, second (DT) as contourplot


i=1;
j=3;
plotdFMcc(name, 1,1,FIM_TS,FIM_DT,i,j,'TS','DT'); 
plotdFM(name, 1,1,FIM_TS,FIM_DT,i,j);



plotdFMAllcc(name, 1,1,FIM_TS,FIM_DT,'TS','DT')


% sensitivity coeffitients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SC1=sensitivities(F{1});


sensitivitiesAll(name, F);

diaginv(name,F);



decomp(name,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

addpath(genpath('lib'))
name='singlegene';

create(name);
addpath(genpath(['models/',name]))

y0(1)=50;
y0(2)=1000;
y0(3)=1;
y0(4)=1;
y0(5)=0;

plottraj(name,10,0.3,5,y0,[1, 2])
F=Fisher(name,10,0.1,3,y0,[1 2],10,'TS','FALSE');
FIM_TS=F;
F=Fisher(name,10,0.1,3,y0,[1 2],10,'TP','FALSE');
FIM_TP=F;
F=Fisher(name,10,0.1,3,y0,[1 2],10,'DT','FALSE');
FIM_DT=F;

[FIM_All]=Fisher(name,10,0.1,3,y0,[1 2],10,'All','FALSE');


F=Fisher(name,10,0.1,3,y0,[1 2],10,'TS','TRUE');
FIM_TSlog=F;
F=Fisher(name,10,0.1,3,y0,[1 2],10,'TP','TRUE');
FIM_TPlog=F;
F=Fisher(name,10,0.1,3,y0,[1 2],10,'DT','TRUE');
FIM_DTlog=F;

[FIM_Alllog]=Fisher(name,10,0.1,3,y0,[1 2],10,'All','TRUE');



F=Fisher(name,10,0.1,3,y0,[1 2],10,'DT','FALSE');
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

clear all;
close all;
name='EnzStat';
create(name);
addpath(genpath(['models/',name]))

y0(1)=0;
y0(2)=120;
y0(3)=0;
y0(4)=0;

y0(5)=10;
y0(6)=10;
y0(7)=10;
y0(8)=10;

y0(9)=0;
y0(10)=0;
y0(11)=0;

y0(12)=0;
y0(13)=0;
y0(14)=0;

 

plottraj(name,10,1,15,y0,[1 4])



F=Fisher(name,10,0.1,3,y0,[1 4],10,'All','T');
FIM_TS=F{1};
FIM_TP=F{2};
FIM_DT=F{3};



plotdFMAll(name, 1,1,FIM_TS,FIM_DT); % First element plotted (TS) polted as heatmap, second (DT) as contourplot


i=1;
j=3;
plotdFMcc(name, 1,1,FIM_TS,FIM_DT,i,j,'TS','DT'); 




plotdFMAllcc(name, 1,1,FIM_TS,FIM_DT,'TS','DT')


% sensitivity coeffitients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SC1=sensitivities(F{1});


sensitivitiesAll(name, F);

diaginv(name,F);



decomp(name,F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = plottraj(name,N,freq, init_T,y0,obs)

[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
addpath([pwd,'/models/','/',name,'/symbolic/']);
addpath([pwd, '/models','/',name]);




[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
size_sto=size(stoichiometry); 
nvar=size_sto(1);
par=parv;
npar=length(parn);
all_equations=str2func([name,'_all_equations']);
ext_nvar=2*nvar+(nvar-1)*nvar/2;
x = linspace(init_T,init_T+freq*N,N);
tspan=[0,init_T+freq*N];

mysolution=ode15s(all_equations,tspan,y0,[],parv);
traj=deval(mysolution,x);

VAR=sqrt(traj((nvar+1):(2*nvar),:));


figure
ii=1
for i=obs,
 subplot(length(obs),1,ii);
 errorbar(x,traj(i,:)',VAR(i,:)') 
legend((['variable  ', int2str(i)]));
ii=ii+1;
end


suptitle('Model trajectories with their standard devaiations')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









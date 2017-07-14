function R = plot_io(name,N,freq,init_T,y0,obs,par,input_ind,output_time)

addpath(genpath([pwd, '/models','/',name]));

stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
size_sto=size(stoichiometry);
size_par=size(par);
nvar=size_sto(1);
npar=size_par(2);
nbatch=size_par(1);

all_equations=str2func([name,'_all_equations']);
nvar_ext=2*nvar+(nvar-1)*nvar/2;

x = linspace(init_T,init_T+freq*(N-1),N);
tspan=[0,init_T+freq*(N-1)];

traj_mean=zeros(nvar_ext,N,nbatch);
traj_std=zeros(nvar,N,nbatch);

for k=1:nbatch
	mysolution 	= ode15s(all_equations,tspan,y0,[],par(k,:) );
	traj_mean(:,:,k) 	= deval(mysolution,x);
	traj_std(:,:,k) 	= sqrt(traj_mean(((nvar+1):(2*nvar)),:,k));
end

x_input=par(:,input_ind);

figure;set(gcf,'Visible', 'off'); 
for i=1:length(obs),
 subplot(length(obs),1,i);
 plot(x_input,  squeeze(traj_mean(obs(i),output_time,:)) ); 
 legend((['variable  ', int2str(i)]));
end

suptitle('Model trajectories with their standard devaiations')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
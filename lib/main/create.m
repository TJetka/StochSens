function [ output_args ] = create(name)

%% reading the model definition files 
disp('reading the model definition files') 
DefsDir = [pwd, '/models/', name];  %where model definition files are found
[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); %reading in information about parameters
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']); %reading stoichiometry matrix
mkdir([pwd,'/models/','/',name,'/symbolic/']); %creating a folder for results of symbolic computations 
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding paths to the model
addpath([pwd, '/models','/',name]);
ratestxt = fileread([pwd, '/models/','/',name,'/' ,name,'_rates.m']);  %name of the function contiaing reaction rates

size_sto=size(stoichiometry); %dimension of stoichiomety matrix
varn=size_sto(1); % number of variables
rean=size_sto(2); % number of reactions

%%symbolic definition of stoichiometry matrix
syms S;
for i=1:size_sto(1),
    for j=1:size_sto(2),
    S(i,j)=stoichiometry(i,j);
    end
end
 
%%symbolic definition of the state vector
disp('creating definition of the model');
y=sym('y','positive');
% create names of variables

for i=1:varn,   
     y(i) = sym((['y',num2str(i)]),'positive');
     var_st{i}=['y(',num2str(i),')'];
     var_st_sym{i}=['y',num2str(i)];
end

size_stim=max(cell2mat(cellfun(@(x) str2num(x{1}),regexp(regexp(ratestxt,'stimulus\([0-9]+\)','match'),'[0-9]+','match'), 'UniformOutput', false)));

% symbolic definition of the stimulus
stim=sym('stim','positive' );
for i=1:1:size_stim
    stim(i)=sym(['stim',num2str(i)],'real' );
    stm_st{i}=[name,'_stimulus(t,',num2str(i),')'];  %name using brackets
    stm_st_sym{i}=['stim',num2str(i)]; %name without brackets
end

t=sym('t','positive' );

%%symbolic definition of parameter vector
param=sym('param','positive' );

pnum = length(parn);
for i = 1:pnum,
     param(i) = sym(['p',num2str(i)],'positive'); 
     par_st{i}=['p(',num2str(i),')'];
     par_st_sym{i}=['p',num2str(i)];
end
rates = str2func([name,'_rates']);

%% creating  symbolic macroscopic rate equation (MRE)
disp('creating macroscopic rate equation')
if length(regexp(ratestxt,'t[ ]*,[ ]*stimulus','match'))==0
    disp('warning:no stimulus detected');
    symrates=feval(rates,y,param);
elseif length(regexp(ratestxt,'t[ ]*,[ ]*stimulus','match'))>0
    symrates=feval(rates,y,param,t,stim);
else
   error('something went wrong in create.m'); 
end
MRE=S*symrates;

%% creating symbolic Jacobian of MRE
 disp('creating Jacobian of MRE')
 J=jacobian(MRE,y(1:varn));
 J_pom=substitution(J,var_st_sym,var_st);
 if size_stim>0 
     J_pom=substitution(J_pom,stm_st_sym,stm_st);
 end 
savefunction(substitution(J_pom,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_MRE_jacobian.m']);


%% creating symbolic variance equations
disp('creating variance equations')
E=sym('E','positive');
for j=1:rean,
E(j,j)=(sqrt(symrates(j)));        
end
EE=sym('EE','positive');
for j=1:rean,
EE(j,j)=((symrates(j)));        
end
D=sym('D','real');
D=(S*EE)*(S');
Sigma=sym('Sigma','real');

totdim=2*varn+varn*(varn-1)/2;

y(totdim)=0;

% giving names to symbolic variables
k=1;
for i=1:varn,
for j=(i+1):varn,
     y(varn+varn+k) = sym(['y',num2str(varn+varn+k)],'real');
     var_st{varn+varn+k} = ['y(',num2str(varn+varn+k),')'];
     var_st_sym{varn+varn+k} = ['y',num2str(varn+varn+k)];
 k=k+1;   
end
end

for i=1:varn,
     y(varn+i) = sym(['y',num2str(varn+i)],'real');
     var_st{varn+i}=['y(',num2str(varn+i),')'];
     var_st_sym{varn+i}=['y',num2str(varn+i)];
end

% assigning names to symbolic matrix Sigma
for i=1:varn,%diagonal elements
  Sigma(i,i)=y(varn+i);
end

%%% lower diagonal
k=1;
for i=1:varn,
for j=(i+1):varn,
   Sigma(i,j)=y(varn+varn+k);
 k=k+1;   
end
end
%using symmetry
Sigma=Sigma+Sigma';

%%%diagonal
for i=1:varn,
   Sigma(i,i)=y(varn+i);
end

Sigma_dot=J*Sigma+ Sigma*(J')+D;

% rewriting variance equations from a matrix format to a vector format
sym('variances_dot','real');

k=1;
for i=1:varn,
variances_dot(i)=Sigma_dot(i,i);
k=k+1;
end

for i=1:varn,
    for j=(i+1):varn,
     variances_dot(k)=Sigma_dot(i,j);
     k=k+1;   
    end
end

%% create concatenated set of equations Y
all_eq=[MRE; variances_dot'];
all_eq_pom=substitution(all_eq,var_st_sym,var_st);
 if size_stim>0 
    all_eq_pom=substitution(all_eq_pom,stm_st_sym,stm_st);
 end
%writting equation as a function
savefunction(substitution(all_eq_pom,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_all_equations.m']);

disp('creating Jacobians')
%% create variable jacobian of Y
J_all_eq_dvar=jacobian(all_eq,y);
J_all_eq_dvar=substitution(J_all_eq_dvar,var_st_sym,var_st);
 if size_stim>0 
    J_all_eq_dvar=substitution(J_all_eq_dvar,stm_st_sym,stm_st);
 end
%writting equation as a function
savefunction(substitution(J_all_eq_dvar,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_all_equations_jacobian_dvar.m']);

%% create parameter jacobian of Y
J_all_eq_dpar=jacobian(all_eq,param); % variables in raws parameters in colums;each row is differentiated by each parameters subsequently
J_all_eq_dpar=substitution(J_all_eq_dpar,var_st_sym,var_st);
 if size_stim>0 
    J_all_eq_dpar=substitution(J_all_eq_dpar,stm_st_sym,stm_st);
 end
%writting equation as a function
savefunction(substitution(J_all_eq_dpar,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_all_equations_jacobian_dpar.m']);


%% variable jacobians  of variable jacobian of MRE
mkdir([pwd, '/models/','/',name,'/symbolic/J_JMRE_dvar'])
for i=1:varn
    J_JMRE_dvar(:,i,:)=jacobian(J(:,i),y(1:varn));
end

for i=1:varn
    %writting equation as a function
    J_JMRE_dvar(:,:,i)=substitution(J_JMRE_dvar(:,:,i),var_st_sym,var_st);
     if size_stim>0 
        J_JMRE_dvar(:,:,i)=substitution( J_JMRE_dvar(:,:,i),stm_st_sym,stm_st);
     end
    savefunction(substitution(J_JMRE_dvar(:,:,i),par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/J_JMRE_dvar/',name,'_jacobianMRE_jacobian_dvar_',num2str(i),'.m']);
end

%% parameter jacobians  of variable jacobian of MRE
mkdir([pwd, '/models/','/',name,'/symbolic/J_JMRE_dpar'])
for i=1:varn
    J_JMRE_dpar(:,i,:)=jacobian(J(:,i),param);
end

for i=1:pnum
    J_JMRE_dpar(:,:,i)=substitution(J_JMRE_dpar(:,:,i),var_st_sym,var_st);
     if size_stim>0 
        J_JMRE_dpar(:,:,i)=substitution(J_JMRE_dpar(:,:,i),stm_st_sym,stm_st);
     end
    %writting equation as a function
    savefunction(substitution(J_JMRE_dpar(:,:,i),par_st_sym,par_st), [pwd, '/models/','/',name,'/symbolic/J_JMRE_dpar/',name,'_jacobianMRE_jacobian_dpar_',num2str(i),'.m']);
end

disp('files created')
end
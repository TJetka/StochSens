function [FisherTS, FisherTP, FisherDT ]= FisherFunctionLogTS(name,N,freq, init_T,y0,obs,merr)
% function calculates FIM for TS data for standard parameters

%%path to a folder where model files are stored
addpath([pwd, '/models','/',name,'/symbolic/' ]);

%% read in parameters of a  model
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
%% read in stoichiometry of a  model    
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);
npar=length(parn);
nvar=size(stoichiometry);
nvar=nvar(1);
ext_nvar=2*nvar+(nvar-1)*nvar/2;

%% defining references to model dependent functions created by create('modelname')  
global all_equations, 
all_equations=str2func([name,'_all_equations']); % mean and variance equation
all_equations_jacobian_dpar=str2func([name,'_all_equations_jacobian_dpar']); % parameter Jacobian of mean and variance equation
all_equations_jacobian_dvar=str2func([name,'_all_equations_jacobian_dvar']); % variables Jacobian of mean and variance equation
MRE_jacobian=str2func([name,'_MRE_jacobian']); % variable Jacobian of mean equation

for i=1:nvar
  jacobianMRE_jacobian_dvar{i}=str2func([name,'_jacobianMRE_jacobian_dvar_',num2str(i)]); % variable  Jacobian of Jacobian of mean equation
end

for i=1:npar
  jacobianMRE_jacobian_dpar{i}=str2func([name,'_jacobianMRE_jacobian_dpar_',num2str(i)]); % parameter Jacobian of Jacobian of mean equation
end

%% calculate trajetories i.e. Y from the paper - concatenation of means and variances 
disp('Calculating trajectories '); 
tspan=[0,init_T+freq*(N-1)]; %time range
x = linspace(init_T,init_T+freq*(N-1),N); %time
mysolution=ode15s(all_equations,tspan,y0,[],par); %solving Y equation
traj=deval(mysolution,x); % evaluating trajectory

%% calulate derivatives of Y with respect to each parameter
disp('Calculating derivatives')
for i=1:npar, % for each parameter
  traj_derivative(i)= ode15s(@der_dpark,tspan, zeros(ext_nvar,1),[],mysolution, par,i,all_equations_jacobian_dpar, all_equations_jacobian_dvar); %solving appropriate ODE
  traj_derivative_val(:,:,i)=par(i)*deval(traj_derivative(i),x); %evaluating solution
end

%%  Exctarct covariance matrix at each observation time point 
Sigma_t=zeros(nvar,nvar,N); % matrix containing variable covariances for all nvar variables

% extracting values from ODE solution using the fact that Sigma is symmetric
tempMatrix=zeros(nvar,nvar);
for k=1:N,
  tempMatrix(tril(true(nvar),-1))=traj((2*nvar+1):end,k);
  Sigma_t(:,:,k)=(tempMatrix)+(tempMatrix)' + diag(traj( (nvar+1):(2*nvar),k));
end

% for k=1:N,
%   l=2*nvar+1;    
%   for i=1:nvar,
%     for j=(i+1):nvar,
%       Sigma_t(i,j,k)=traj(l,k);  
%       l=l+1;
%     end
%   end
% end

% for k=1:N,
%   l=nvar+1;    
%   Sigma_t(:,:,k)=Sigma_t(:,:,k)+Sigma_t(:,:,k)';
%   for i=1:nvar,
%     Sigma_t(i,i,k)=traj(l,k);
%     l=l+1;
%   end
% end

%% calculate derivative of a covariance matrix at each observation Sigma_der_t
Sigma_der_t=zeros(nvar,nvar,N,npar);

tempMatrix=zeros(nvar,nvar);
for r=1:npar, %derivative index
  for k=1:N,  %time-point index
    tempMatrix(tril(true(nvar),-1))=traj_derivative_val((2*nvar+1):end,k,r);
    Sigma_der_t(:,:,k,r)=(tempMatrix)+(tempMatrix)' + diag(traj_derivative_val( (nvar+1):(2*nvar),k,r));
  end
end

% for r=1:npar, %derivative index
%   for k=1:N,  %time-point index
%     l=2*nvar+1;    
%     for i=1:nvar, %i,j - covariance matrix index
%       for j=(i+1):nvar,
%         Sigma_der_t(i,j,k,r)=traj_derivative_val(l,k,r);   %extracting from ODE solution  
%       end
%     end
%   end

%   for k=1:N,%time-point index
%     l=nvar+1;    
%     Sigma_der_t(:,:,k,r)=Sigma_der_t(:,:,k,r)+Sigma_der_t(:,:,k,r)';
%     for i=1:nvar, %diagonal index
%       Sigma_der_t(i,i,k,r)=traj_derivative_val(l,k,r);
%       l=l+1;
%     end
%   end
% end

%% Fundamental matrices
disp('Calculating fundamental matrices');
Fund=zeros(nvar,nvar,N,N); %defining the 2 dimensional list of matrices 
Fund_traj{nvar}=[]; 
for i=1:(N-1), % for time point i=1,...,N where we start from
  tspanF=[x(i),x(N)]; %tme range 
  for j=1:nvar,% for each canonical vector
    ei=zeros(nvar,1);
    ei(j)=1;
    Fund_sol=ode15s(@FundRhs,tspanF,ei,[],par, mysolution,nvar, MRE_jacobian); %solving appropriate ODE
    Fund_traj{j}=deval(Fund_sol,x((i+1):N)); %evaluating
  end
  for j=1:nvar
    POM=Fund_traj{j};
    kk=1;
    for ii=(i+1):N, % where we go to 
      Fund(:,j,i,ii)=POM(:,kk);% i where we start from ii where we got to  
      kk=kk+1;
    end     
  end
end

%% building large covariance matrix 
SigmaBlocks=zeros(nvar,nvar,N,N);

for i=1:(N),
  SigmaBlocks(:,:,i,i)=Sigma_t(:,:,i); %assigning diagonal elements
  for ii=(i+1):N,
    SigmaBlocks(:,:,i,ii)= (Fund(:,:,i,ii)*SigmaBlocks(:,:,i,i))'; % filling rows 
  end         
end

%% Fundamental matrices derivatives
dFund=zeros(nvar,nvar,N,N,npar); %defining 3 dimensional list of matrices
Fund_traj_der{nvar}=[];

for r=1:npar,
  for i=1:(N-1),
    tspanF=[x(i),x(N)];
    
    for j=1:nvar, % for each canonical vector
      ei=zeros(nvar,1);
      zz=zeros(nvar,1);
      ei(j)=1;  
      Fund_sol=ode15s(@FundRhs,tspanF,ei,[],par, mysolution,nvar,MRE_jacobian); %Solving for fundamental matrix
      Fund_sol_der=ode15s(@derJ_dpark,tspanF,zz,[], mysolution,traj_derivative, Fund_sol, par,nvar,r, MRE_jacobian, jacobianMRE_jacobian_dpar,  jacobianMRE_jacobian_dvar); %Solving for derivative
      %log derivative
      Fund_traj_der{j}=par(r)*deval(Fund_sol_der,x((i+1):N)); %evaluating  
    end
    
    for j=1:nvar,
      POM=Fund_traj_der{j};
      kk=1;
      for ii=(i+1):N,
        dFund(:,j,i,ii,r)=POM(:,kk); % populating the 3 dimensional list of derivatives of fundamental matrices
        kk=kk+1;
      end    
    end

  end
end

dSigmaBlocks=zeros(nvar,nvar,N,N,npar); % defining  3 dimensional list of derivatives of variances

for r=1:npar, % for each parameter
  for i=1:(N), % for each time point
    dSigmaBlocks(:,:,i,i,r)=Sigma_der_t(:,:,i,r); % populating diagonal elements
    for ii=(i+1):N, %populating rows
      dSigmaBlocks(:,:,i,ii,r)= (Fund(:,:,i,ii)*dSigmaBlocks(:,:,i,i,r)+dFund(:,:,i,ii,r)*SigmaBlocks(:,:,i,i))';
    end    
  end
end

%% large covariance matrix and its derivatives
Sigma=zeros(nvar*N,nvar*N);  % covariance matrix of observed variables
first=1;
last=nvar;
    
first_c=nvar+1;
last_c=first_c+nvar-1;
    
% popuulating upper diagonal covariance matrix from list SigmaBlocks
for i=1:(N-1),
  %Sigma(first:last,first:last)=Sigma_t(:,:,i);
  first_c=last+1;
  last_c=first_c+nvar-1;
  for ii=(i+1):N,       
    Sigma(first:last,first_c:last_c)=SigmaBlocks(:,:,i,ii);
    first_c=last_c+1;
    last_c=first_c+nvar-1;
  end    
  first=last+1;
  last=first+nvar-1;     
end

% populating the rest of the covariance matrix 
Sigma=Sigma+Sigma';
first=1;
last=nvar;

for i=1:N,
  Sigma(first:last,first:last)=SigmaBlocks(:,:,i,i)+diag(merr*ones(1,nvar));
  first=last+1;
  last=first+nvar-1;   
end

%% derivative of large covariance matrix
dSigma=zeros(nvar*N,nvar*N,npar);

%results kept in dSigmaBlocks will be put into the large 2d mattrix dSigma
for r=1:npar, %for each parameter
  first=1;
  last=nvar;
  first_c=nvar+1;
  last_c=first_c+nvar-1;
  for i=1:(N-1), % for each time point except the last one  
    %  Sigma(first:last,first:last)=Sigma_t(:,:,i);
    first_c=last+1;
    last_c=first_c+nvar-1;
    for ii=(i+1):N,
      dSigma(first:last,first_c:last_c,r)=dSigmaBlocks(:,:,i,ii,r);
      first_c=last_c+1;
      last_c=first_c+nvar-1;
    end    
    first=last+1;
    last=first+nvar-1;     
  end
end
    
for r=1:npar,
  dSigma(:,:,r)=dSigma(:,:,r)+dSigma(:,:,r)';   
end

for r=1:npar,
  first=1;
  last=nvar;
  for i=1:N,
    dSigma(first:last,first:last,r)=dSigmaBlocks(:,:,i,i,r);
    first=last+1;
    last=first+nvar-1; 
  end
end

%% Concatenetions
% concatenation of calculated derivatives into format convenient to
% evaluate FIM
conc_traj_deriv=zeros(N*nvar,npar);

for r=1:npar,
  first=1;
  last=nvar;
  for k=1:N,    
    conc_traj_deriv(first:last,r)=traj_derivative_val(1:nvar,k,r);
    first=last+1;
    last=first+nvar-1; 
  end
end

Sigma=ConcMat(Sigma,obs,nvar,N);
conc_traj_deriv=ConcVec(conc_traj_deriv,obs,nvar,N);
 
lobs=length(obs);

dSigma_pom=dSigma;
dSigma=zeros(lobs*N,lobs*N,npar);

for i=1:npar,  
  dSigma(:,:,i)=ConcMat(dSigma_pom(:,:,i),obs,nvar,N);
end

S=Sigma;
[L U P]=lu(S);
Fisher=zeros(npar,npar);
Var_term=zeros(npar,npar);
Mean_term=zeros(npar,npar);

% Using LU factorization
for i=1:npar,
  Ci=conc_traj_deriv(1:(N*lobs),i);
  Di=dSigma(1:(N*lobs),1:(N*lobs),i);
  for j=(i+1):npar,
    Cj=conc_traj_deriv(1:(N*lobs),j);
    Dj=dSigma(1:(N*lobs),1:(N*lobs),j);
    Mean_term(i,j)=Ci'*(U\(L\(P*Cj)));
    Var_term(i,j)=0.5*trace((U\(L\(P*Di)))*(U\(L\(P*Dj))));
    Fisher(i,j)= Mean_term(i,j)+Var_term(i,j);  
  end
end

Fisher=Fisher+Fisher';

for i=1:npar,
  Ci=conc_traj_deriv(1:(N*lobs),i);
  Di=dSigma(1:(N*lobs),1:(N*lobs),i);
  Mean_term(i,i)=Ci'*(U\(L\(P*Ci)));
  Var_term(i,i)=0.5*trace((U\(L\(P*Di)))*(U\(L\(P*Di))));
  Fisher(i,i)= Mean_term(i,i)+Var_term(i,i);
end
FisherTS=Fisher;

clear global all_equations  all_equations_jacobian_dpar all_equations_jacobian_dvar  MRE_jacobian jacobianMRE_jacobian_dpar jacobianMRE_jacobian_dvar
end
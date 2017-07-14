function OUT  = ConcMat(Sigma,obs,nvar,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sizig=size(Sigma);

nobs=length(obs);


OUT=zeros(N*nobs, N*nobs);


begO_r=1;
endO_r=nobs;



begS_r=1;
endS_r=nvar;


begS_c=1;
endS_c=nvar;


for i=1:N,

    
begS_c=1;
endS_c=nvar;


begO_c=1;
endO_c=nobs;

    
    
for j=1:N,
    
   O_r=begO_r:endO_r; 
   O_c=begO_c:endO_c; 
    
   
   S_r=begS_r:endS_r; 
   S_c=begS_c:endS_c; 
   
   OUT(O_r,O_c)=Sigma(S_r(obs),S_c(obs));
   
   begO_c=endO_c+1;
   endO_c=begO_c+nobs-1;
   
   
   begS_c=endS_c+1;
   endS_c=begS_c+nvar-1;
 
end

  begO_r=endO_r+1;
  endO_r=begO_r+nobs-1;
  begS_r=endS_r+1;
  endS_r=begS_r+nvar-1;
  
end


end


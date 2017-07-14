function [ OUT ] = ConcVec(Vec,obs,nvar,N)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

sizeV=size(Vec);
nobs=length(obs);
OUT=zeros(N*nobs, sizeV(2));

begO_r=1;
endO_r=nobs;
begS_r=1;
endS_r=nvar;


for i=1:N,
    
   O_r=begO_r:endO_r; 
   S_r=begS_r:endS_r; 
   
   OUT(O_r,:)=Vec(S_r(obs),:);
   
   begO_r=endO_r+1;
   endO_r=begO_r+nobs-1;
   begS_r=endS_r+1;
   endS_r=begS_r+nvar-1;
      
end



end
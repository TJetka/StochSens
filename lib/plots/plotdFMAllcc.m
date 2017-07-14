function [ output_args ] = plotdFMAllcc(name,rx,ry,C,D,L1,L2)

np=size(C);
np=np(1);
kk=1;
figure
for j=1:np,
for i=1:np,
subplot(np,np,kk);
plotdFMcc(name, rx,ry,C,D,j,i,L1,L2)
kk=kk+1;
end
end

end


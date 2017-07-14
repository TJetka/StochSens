function [ output_args ] = plotdFMAll(name,rx,ry,C,D)

np=size(C);
np=np(1);
kk=1;
figure
for j=1:np,
for i=1:np,
subplot(np,np,kk);
plotdFM(name, rx,ry,C,D,j,i)
kk=kk+1;
end
end

end


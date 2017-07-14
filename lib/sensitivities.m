function SC = sensitivities(FM)

sz=size(FM);
[B C]=eig(FM);

Sens=zeros(1,sz(1));

NT=C*(B');

for i=1:sz,
        Sens(i)=norm(NT(:,i));
end
    
Sens=Sens/(sum(Sens));
SC=Sens;
end


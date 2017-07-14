function SC = decomp(name,F)
[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');

sz=size(F{1});
sz=sz(1);

[A B C]=svd(F{1});

C1=C';
SS=zeros(sz,sz);


[A B C]=svd(F{2});
C2=C';

SSTP=zeros(sz,sz);

[A B C]=svd(F{3});

C3=C';
SSDT=zeros(sz,sz);

for i=1:sz,
for j=1:sz,
    SS(i,j)=C1(i,j)^2;
    SSTP(i,j)=C2(i,j)^2;
    SSDT(i,j)=C3(i,j)^2;
end
end

decomp1=figure;
bar3(SS');
set(gca, 'YTick', 1:sz, 'YTickLabel', parn);
ylabel('parameters')
xlabel('eigenvalues')
title('TS');




decomp2=figure;
bar3(SSTP');
ylabel('parameters')
xlabel('eigenvalues')
title('TP');
set(gca, 'YTick', 1:sz, 'YTickLabel', parn);


decomp3=figure;
bar3(SSDT');
%set(gca, 'YTick', 1:sz, 'YTickLabel', parn);
ylabel('parameters')
xlabel('eigenvalues')
title('DT');
set(gca, 'YTick', 1:sz, 'YTickLabel', parn);


end


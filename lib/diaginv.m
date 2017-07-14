function SC = diaginv(name,F)
[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');


FMinv=inv(F{1});
FMTPinv=inv(F{2});


FMinvDiag=diag(FMinv);
FMTPinvDiag=diag(FMTPinv);


BART=figure;
bar([FMinvDiag';FMTPinvDiag']')

set(gca, 'XTick', 1:7, 'XTickLabel', parn);
  legend('TS','TP');
xlabel('parameters');
ylabel('inverse diagonal elements') ; 
  title('inverse  diagonals');
end


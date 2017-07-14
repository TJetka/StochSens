function SC = sensitivitiesAll(name,F)
[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');
SC1=sensitivities(F{1});
SC2=sensitivities(F{2});
SC3=sensitivities(F{3});



BART=figure;
bar([SC1;SC2;SC3]')
set(gca, 'XTick', 1:length(parv), 'XTickLabel', parn);

legend('TS','TP','DT');
xlabel('parameters');
ylabel('sensitivities') ; 
  title('Overall sensitivities');
end


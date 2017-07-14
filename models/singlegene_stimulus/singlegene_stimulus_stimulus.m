function result = JAKSTAT_stimulus(t,i)
    stim={};
    stim{1}=10*(t>1&t<2);

    result=stim{i};
end


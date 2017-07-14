function result = JAKSTATtest_stimulus(t,i)
    result=[];
    stim{1}=0.01*(5*t+1+sin(t));
    stim{2}=20;
    stim{3}=1*((t>4)&(t<8));
    result=stim{i};
end


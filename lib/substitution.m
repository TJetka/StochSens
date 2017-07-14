function [dif, s1, t1] = substitution(dif1, vp, vps);

t1={};
s1={};

if dif1 == 0
    dif = dif1;
    return
end

s = findsym(dif1);

rem = s;
 

i = 0;
while length(rem)
    [t, rem] = strtok(rem,',');
    t = strtok(t);  
    for j=1:length(vp)            
        if strcmp(t,vp{j})
            i=i+1;
            t1{i} = t;
            s1{i}=vps{j};
            break;                   
        end
    end
end

if length(t1) 
    dif = subs(dif1,t1,s1); 
else
    dif = dif1;
end

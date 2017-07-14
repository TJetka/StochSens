function  F = Fisher(name,N,freq, init_T,y0,obs,merr, type, LOG)

%this is an interface function that allows a unified access to specific
%function for FIM compuatation

switch type
    
    case('TP')
       
        if(strcmp(LOG,'TRUE'))
        % TP data for log parameters
        F=FisherFunctionLogTP(name,N,freq,init_T,y0,obs,merr);
        else
        % TP data for standard parameters
        F=FisherFunctionTP(name,N,freq, init_T,y0,obs,merr);
        end
        
    case('TS')   
        if(strcmp(LOG,'TRUE'))
        % TS data for log parameters
        F=FisherFunctionLogTS(name,N,freq,init_T,y0,obs,merr);
        else    
       % TS data for standard parameters       
        F=FisherFunctionTS(name,N,freq, init_T,y0,obs,merr);
        end
        
    case('DT')
        
        if(strcmp(LOG,'TRUE'))
        % DT data for log parameters
        F=FisherFunctionLogDT(name,N,freq,init_T,y0,obs,merr);
        else
        % DT data for standard parameters
        F=FisherFunctionDT(name,N,freq, init_T,y0,obs,merr);
        end
    
    case('All')
   
        if(strcmp(LOG,'TRUE'))          
        % TS, TP and DT  data for log parameters
        [F{1} F{2} F{3}]=FisherFunctionLogAll(name,N,freq, init_T,y0,obs,merr);        
        else
        % TS, TP and DT  data for standard parameters
        [F{1} F{2} F{3}]=FisherFunctionAll(name,N,freq, init_T,y0,obs,merr);
            
        end
end

end


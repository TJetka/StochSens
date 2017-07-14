function R = FundRhs(t,y,p,solution,nvar, MRE_jacobian)

%MRE_jacobian

traj=deval(solution,t);

R=feval(MRE_jacobian,t,traj(1:nvar),p)*y;

end


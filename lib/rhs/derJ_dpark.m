function rhs = derJ_dpark(t,y0,mysolution,traj_derivative,Fund_sol, p,nvar,k,MRE_jacobian,  jacobianMRE_jacobian_dpar,  jacobianMRE_jacobian_dvar)
	%dydt = jac * Y + bs'; %this should convert the rhs matrix into a vector
	% global MRE_jacobian  jacobianMRE_jacobian_dpar  jacobianMRE_jacobian_dvar
	traj=deval(mysolution,t);
	traj_derpar=deval(traj_derivative(k),t);
	traj=traj(1:nvar);
	traj_derpar=traj_derpar(1:nvar);

	Fund=deval(Fund_sol,t);
	J=feval(MRE_jacobian,t,traj,p);
	bs= feval(jacobianMRE_jacobian_dpar{k},t,traj,p);
	
	Kaux=zeros(nvar,nvar);
	for i=1:nvar
		Kaux=Kaux+feval(jacobianMRE_jacobian_dvar{i},t,traj,p)*traj_derpar(i);
	end

	%bt*traj_derpar
	K=(bs + Kaux )*Fund;
	rhs=J*y0+K;
end
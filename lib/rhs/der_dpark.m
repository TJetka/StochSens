function rhs = der_dpark(t,y0,odesol,p,k,all_equations_jacobian_dpar,  all_equations_jacobian_dvar)
	%dydt = jac * Y + bs'; %this should convert the rhs matrix into a vector
	%global all_equations_jacobian_dpar  all_equations_jacobian_dvar

	traj=deval(odesol,t);
	bs = feval(all_equations_jacobian_dpar,t,traj,p);        % evaluates dF/dk
	bs = bs(:,k);               %new
	J=feval(all_equations_jacobian_dvar,t,traj,p);
	rhs=J*(y0)+bs;

end
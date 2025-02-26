function out=linear_solve(inlet_nodes,outlet_nodes,alpha,beta,R,a,qin)

[p,q] = linear_analysis_pressure(inlet_nodes,outlet_nodes,alpha,beta,R,a);

out = qin - q(2);

end

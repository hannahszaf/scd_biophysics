
function [p,q, A]=linear_analysis_pressure(inlet_nodes,outlet_nodes,alpha,beta,R,a)
%BCs
N=max(max(alpha,beta));
%Initialise A
A=eye(N);
%Apply conservation of mass over all nodes
for i=1:N
    %ignore outlets
    if any(outlet_nodes==i)==1
        continue
    end
    %fix PRESSURE at inlets
    if any(inlet_nodes==i)==1
        A(i,i)=1;
        %A(i,beta(find(alpha==i)))=-1;  this is ignored unless fixing Qin
        continue
    end
    %store connection nodes
    c1=find(alpha==i);
    c2=find(beta==i);
    %sum of inverse of connecting resistances
    K=sum(R(c1).^(-1))+sum(R(c2).^(-1));
    %Note identities of nonzero entries in matrix are nodes that LINK to
    %node i
    for j=c1
        A(i,beta(j))=-R(j)^(-1)/K;
    end
    for j=c2
        A(i,alpha(j))=-R(j)^(-1)/K;
    end
end

%Solve for pressuure
p=A\a;
%Generate flow rates
for i=1:size(beta,2)
    q(i)=(p(alpha(i))-p(beta(i)))/R(i);
end
end
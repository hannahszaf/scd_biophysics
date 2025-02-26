function [Rbulk,Rslip,vav,vslip,P, Q, Qtest, res, rexp, Pnew, A, Res2, Pressure2, Rbulk_alg, Rslip_alg, R1, R3, Q_21, Pressure_bypass,Rslip_conQ, Rbulk_conQ, Q_21_sum, Q_byp, Q_byp_lin, Rbulk_oxy, Rslip_oxy, P2_oxy, Qbulk, Qslip, g] = calculate_properties_fit_final(cslip,asym,blunt,vmax,alpha,beta,inlet_nodes,outlet_nodes,Pin,H,D,L,Ox,i,j)
%vslip = average slip velocity in experimental channel (microns/s)
%vav = average velocity in experimental channel (microns/s)
%alpha, beta define network connections
%inlet_nodes, outlet_nodes define inlet and outlet nodes
%Pin = pressure input (Pa) [note 10342.1 Pa = 1.5 psi]
%Rox = resistances of channels at full oxygenation (Pa s / microns^3)
%Rad = (effective) radius of experimental channel (microns)
%L = length of experimental channel (microns) 
%Ox = vector of parameters (cslip, asym, blunt, vmax) for fully oxygenated flow

[Rox, Q_21, Q_21_sum, Rbulk_oxy, Rslip_oxy, P2_oxy, g] = calculate_properties_ox(Ox,alpha,beta,inlet_nodes,outlet_nodes,Pin,H,D,L); %note relative resistances of inflow and exp/bypass channels needed (at full oxygenation); this is imposed inside the function
res = Rox;
Rbulk_oxy = Rbulk_oxy ;
Rslip_oxy = Rslip_oxy ;
Q_21 = Q_21;
Q_21_sum = Q_21_sum;
% v = @(r) vmax * (1 - (1-cslip)*abs(r-asym).^blunt);
% %calculate vslip
% if asym>0
%     vslip = v(-1);
%     points = [-1:0.00001:asym];
%     points2 = (points - asym)/(1+asym);
% else
%     vslip = v(1);
%     points = [asym:0.00001:1];
%     points2 = (points - asym)/(1-asym);
% end
% vav = 2 * trapz(points2,v(points).*abs(points2));
% vav = v_avg;
% vslip = v_slip;

vslip = cslip*vmax;
vav = blunt;

%Calculate slip flow rate and bulk flow rate from average flow speeds
Qbulk = (vav-vslip)*D*H; %microns^3/s
Qslip = vslip * D*H; %microns^3/s

%Find resistance of experimental channel that gives measured flow rate, and
%corresponding pressure [note this step is only necessary if the inflow
%rate isn't measured at every oxygen tension]
a = [Pin; 0; 0; 0]; %Pressure boundary condition
fun = @(x) linear_solve(inlet_nodes,outlet_nodes,alpha,beta,[Rox(1),x,Rox(3)],a,Qbulk+Qslip);
x0 = [0.001 0.1]; %initial guess of exp channel resistance (in Pa s / microns^3)
Rexp = fzero(fun,x0); %(effective) resistance of exp channel at this oxygen tension (in Pa s / microns^3)
[p,q, A] = linear_analysis_pressure(inlet_nodes,outlet_nodes,alpha,beta,[Rox(1),Rexp,Rox(3)],a); %find pressure at calculated resistance (Pa)
rexp = Rexp;
%solve for R2 and P2 analyitically

R1 = Rox(1);
R3 = Rox(3);
qq = Qbulk+Qslip;

Res2 = (Pin / (qq) - Rox(1))/((Rox(1)/Rox(3) + 1));
Pressure2 = Res2 * (qq);
Q_byp_lin = Pressure2 / R3;
Pressure_bypass = (Q_21(1)-(qq))*R3;
Rslip_conQ = Pressure_bypass / Qslip;
Rbulk_conQ = Pressure_bypass / Qbulk;

Rbulk_alg = Pressure2 / Qbulk;
Rslip_alg = Pressure2 /Qslip;
%Find Rbulk and Rslip from Qbulk and Qslip
Rbulk = p(2)/Qbulk;
Rslip = p(2)/Qslip;
P = p(2);
Pnew = p(:,:);
Q = q(:,:);
Qtest = Qbulk + Qslip;
Q_byp = (Q_21(1)-(Qbulk+Qslip));
A = A;
g = g;
%if j==3
%r = [-1:0.01:1];
%figure
%plot(r,v(r))
%end
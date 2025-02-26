function [Rox,Q_21, Q_21_sum, Rbulk_oxy, Rslip_oxy, P2_oxy, G] = calculate_properties_ox(Ox,alpha,beta,inlet_nodes,outlet_nodes,Pin,H,D,L)
%vslip = average slip velocity in experimental channel (microns/s)
%vav = average velocity in experimental channel (microns/s)
%alpha, beta define network connections
%inlet_nodes, outlet_nodes define inlet and outlet nodes
%Pin = pressure input (Pa) [note 10342.1 Pa = 1.5 psi]
%Rox = resistances of channels at full oxygenation (Pa s / microns^3)
%Rad = (effective) radius of experimental channel (microns)
%L = length of experimental channel (microns) 
%Ox = vector of parameters (cslip, asym, blunt, vmax) for fully oxygenated flow
cslip = Ox(1);
asym = Ox(2);
blunt = Ox(3);
vmax = Ox(4);

G=((12 * 147278 / (1-0.63*(H/40))) * (1/(H^3 * 40))) / ((12 * L / (1-0.63*(D/H))) * (1/(D^3 * H)));   %Jose's device
%G =((12 * 147278 / (1-0.63*(20/40))) * (1/((20^3) * 40))) / ((12 * 13000 / (1-0.971*0.63*(D/H))) * (1/((D^3) * H))); %relative resistance of inflow channel to exp/bif channels at full oxygenation [formulas from Table 3.1 of Bruus, H. (2007). Theoretical Microfluidics, OUP.]

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

vslip = cslip*vmax;
vav = blunt;

%Calculate slip flow rate and bulk flow rate from average flow speeds
Qbulk = (vav-vslip)*D*H; %microns^3/s
Qslip = vslip * D*H; %microns^3/s

%Find resistance of experimental channel that gives measured flow rate, and
%corresponding pressure [note this step is only necessary if the inflow
%rate isn't measured at every oxygen tension]
a = [Pin; 0; 0; 0]; %Pressure boundary condition
fun = @(x) linear_solve(inlet_nodes,outlet_nodes,alpha,beta,[G*x,x,x],a,Qbulk+Qslip);
x0 = [0.001 0.1]; %initial guess of exp channel resistance (in Pa s / microns^3)
Rexp = fzero(fun,x0); %(effective) resistance of exp channel at this oxygen tension (in Pa s / microns^3)
[p,q] = linear_analysis_pressure(inlet_nodes,outlet_nodes,alpha,beta,[G*Rexp,Rexp,Rexp],a); %find pressure at calculated resistance (Pa)
Rox = [G*Rexp, Rexp, Rexp]; %vheck 2
%Find Rbulk and Rslip from Qbulk and Qslip
Rbulk_oxy = p(2)/Qbulk;
Rslip_oxy = p(2)/Qslip;
P2_oxy = p(2);
Q_21 = q;
Q_21_sum = Qbulk + Qslip; 

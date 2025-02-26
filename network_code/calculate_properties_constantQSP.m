%Code to solve Rslip and Rbulk assuming the flowrate stays constant
clear variables
x = csvread('input_UMN030.csv');
%for k = 1:12
i = 1;

%pressure = 2.206*0.8; %CHC041
%pressure = 2.06*0.8; %CHC038
p%ressure = 1.995*0.8; %CHC054
%pressure = 1.995*0.8; %CHC069
%pressure = 2.152*0.8; %CHC072
%pressure = 1.7; %CHC029
%pressure = 1.6; %UMN020
pressure = 2.05*0.8; %UMN030
%x = x(i,:);
%
%x = [0.602600736823780,0.0297360252794870,1.83649699244596,820.969182664471];

%Define network topology
alpha = [1 2 2]; %upstream node IDs
beta = [2 3 4]; %downstream node IDs
inlet_nodes = 1;
outlet_nodes = [3 4];
H = 20.0;
D = 20.0;

%Define parameters
L = 13000; %microns. Length of exp channel
Pin = (pressure) .* 6894.76; %Pressure input in Pa (= 1.5 psi)

%ox_tensions = [0 2 4 6 8 10];
%ox_tensions = [0 2 4 6 10]; %UMN005
%ox_tensions = [0 2 3 4 5 7 12]; %UMN031
%ox_tensions = [0 2 4 6]; %UMN029

%for j = 1:size(ox_tensions,2) %original code
for j = 1:size(x,1)
        %relevant section of dataset
        %y = x(2,(j-1)*4+1:(j-1)*4+4); %origianl code
        %z = [0.735277278 0.039250829 1.352506927 828.6507258]; %signal processing code 21%

       %k = 8; %use for input_fit
       % k = size(ox_tensions,2);
        %z = x(1,(j-1)*4+1:(j-1)*4+4);
        y = x(j,1:7);

        %sweep through samples
        for i = 1 : 1 %size(x,1)
            if isnan(y(i,1))== 0
                %[K(i,j),lambda(i,j),n(i,j),b(i,j)] = calculate_properties(y(i,1),y(i,2),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                %[K(i,j),lambda(i,j),n(i,j),b(i,j),B(i,j)] = calculate_properties_fit(y(i,1),y(i,2),y(i,3),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                %[K(i,j),lambda(i,j),n(i,j),B(i,j)] = calculate_properties_fit_simple(y(i,1),y(i,2),y(i,3),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                [Rbulk(i,j),Rslip(i,j),vav(i,j),vslip(i,j),P(i,j), Q, Qtest(i,j), res, rexp(i,j), Pnew, A, Res2(i,j), Pressure2(i,j), Rbulk_alg(i,j) , Rslip_alg(i,j), R1(i,j), R3(i,j), Q_21, Pressure_bypass(i,j), Rslip_conQ(i,j), Rbulk_conQ(i,j),Q_21_sum, Q_byp(i,j), Q_byp_lin(i,j), Rbulk_oxy(i,j), Rslip_oxy(i,j), P2_oxy(i,j), qbulk(j,i), qslip(j,i), G(i,j)] = calculate_properties_fit_final(y(i,1),y(i,2),y(i,3),y(i,4),y(i,6),y(i,7),alpha,beta,inlet_nodes,outlet_nodes,Pin(i,1),H,D,L,z(i,:),i,j);
            end
        end
end



Reff = 1./((1./Rslip)+(1./Rbulk));
Reff = Reff';
Rslip_conQ = Rslip_conQ'; 
Rbulk_conQ = Rbulk_conQ';
Reff_conQ = 1./((1./Rslip_conQ)+(1./Rbulk_conQ));

Rbulk_oxy = Rbulk_oxy';
Rslip_oxy = Rslip_oxy';
R3 = R3';
P2_oxy = P2_oxy';
Pressure_bypass = Pressure_bypass';
Rfinal = horzcat(P2_oxy,Pressure_bypass,Reff_conQ ,Rbulk_conQ,Rslip_conQ,R3,Rbulk_oxy,Rslip_oxy);

Rbulk = Rbulk';
Rslip = Rslip';


%save('CHC029_Pbypass.mat','Pressure_bypass','P2_oxy');


%save('CHC041_resistance.mat','Rslip','Rbulk','Rslip_oxy','Rbulk_oxy');









%Q stays constant -> Does R change or delP? I think delP...

%delP_bypass = (Q_total - Q_exp)*(R_bypass);

%R_exp = delP_bypass / Q_exp;
%R_slip = delP_bypass / Q_slip;
%R_visc = delP_bypass / Q_bulk;

%need to calculate Q at 21% before each oxygen step 
% we know Qtot, R1, R2, R3, delP at 21%
%0% oxygen we know Qtot, R1, R3, need to find delP and R2
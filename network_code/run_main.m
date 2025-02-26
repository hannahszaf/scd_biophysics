%Code to solve Rslip and Rbulk assuming the flowrate stays constant
%%
clear variables
id = 'input_patient3_new';
id2 = 'pressure_patient3_new';
file = sprintf('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/%s.mat',id);
file2 = sprintf('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/%s.mat',id2);
data = load(file);
data=data.output_params;
data2 = load(file2);
pressure = data2;
pressure = data2.pressure;

date = sprintf('NA');


for k = 1:11
x = squeeze(data(:,:,k)); 
pres = mean(pressure);

res_final(:,:) = calculate_properties_constantQ(x,id,date,pres);
% a = res_final(:,1:3);
% b = res_final(:,4:6);

%res_b(:,k) = res_final(:,1);
res_oxy(:,1,k) = res_final(:,4);
res_oxy(:,2,k) = res_final(:,5);
res_oxy(:,3,k) = res_final(:,6);
% c=a./b;
% res_norm_b(:,k) = c(:,1);
% res_norm_s(:,k) = c(:,2);
% res_norm_e(:,k) = c(:,3);
end
res_oxymean = squeeze(mean(res_oxy,3));

for k = 1:11
x = squeeze(data(:,:,k)); 
pres = mean(pressure);

res_final(:,:) = calculate_properties_constantQ(x,id,date,pres);
% a = res_final(:,1:3);
% b = res_final(:,4:6);

%res_b(:,k) = res_final(:,1);
res_deoxy(:,:) = res_final(:,1:3);

 c=res_deoxy./res_oxymean;
 res_norm_b(:,k) = c(:,1);
 res_norm_s(:,k) = c(:,2);
 res_norm_e(:,k) = c(:,3);
end


res_mean(:,1) = mean(res_norm_b,2); %viscous resistance 
res_mean(:,2) = std(res_norm_b,0,2)./sqrt(k); 
res_mean(:,3) = mean(res_norm_s,2); %friction resistance
res_mean(:,4) = std(res_norm_s,0,2)./sqrt(k); 
res_mean(:,5) = mean(res_norm_e,2); %effective resistance
res_mean(:,6) = std(res_norm_e,0,2)./sqrt(k); 

function res_final = calculate_properties_constantQ(x,id,date,pres)
i = 1;
pressure = pres*0.8;

%pressure = 2.0598*0.8; %input_UMN030
%pressure = 2.0323*0.8; %CHC045
%pressure = 1.995*0.8; %CHC030
%pressure = 2.07*0.8; %CHC058?
%pressure = 2.151*0.8; %UMN027
%pressure = 1.996*0.8; %CHC043
%pressure = 2.0*0.8; %
%pressure = 2.197*0.8; %UMN031
%pressure = 1.5*0.8; %input_MGH2037_TR30
%pressure = 1.8*0.8; %input_MGH2037_TR50
%pressure = 1.92*0.8; %input_MGH2037_TR100
%pressure = 1.75*0.8; %input_CHC004_TR100
%pressure = 1.995*0.8; %input_CHC004_TR30
%pressure = 2.014*0.8; %input_CHC004_TR10
%pressure = 2.6*0.8; %input_CHC058_TR100
%pressure = 1.95*0.8; %input_CHC058_TR10
%pressure = 1.8*0.8; %input_CHC058_TR30
%pressure = 1.9*0.8; %input_CHC058_TR100
%pressure = x(3,1); %maybe put pressure in input data?

%Define network topology
alpha = [1 2 2]; %upstream node IDs
beta = [2 3 4]; %downstream node IDs
inlet_nodes = 1;
outlet_nodes = [3 4];
H = 20.0;
D = 20.0;

%Define parameters
%L = 13000; %microns. Length of exp channel Jose
L = 22000; %microns. Length of exp channel Hannah
Pin = (pressure) .* 6894.76; %Pressure input in Pa (= 1.5 psi)
n=size(x);
ox_tensions = n(2)/4;%[0 2 6];
%ox_tensions = [0 2 4 6 10]; %UMN005
%ox_tensions = [0 2 3 4 5 7 12]; %UMN031
%ox_tensions = [0 2 4 6]; %UMN029

for j = 1:ox_tensions
        %relevant section of dataset
        y = x(2,(j-1)*4+1:(j-1)*4+4);
       %k = 8; %use for input_fit
       % k = size(ox_tensions,2);
        z = x(1,(j-1)*4+1:(j-1)*4+4);
        
        %sweep through samples
        for i = 1 : 1 %size(x,1)
            if isnan(y(i,1))== 0
                %[K(i,j),lambda(i,j),n(i,j),b(i,j)] = calculate_properties(y(i,1),y(i,2),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                %[K(i,j),lambda(i,j),n(i,j),b(i,j),B(i,j)] = calculate_properties_fit(y(i,1),y(i,2),y(i,3),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                %[K(i,j),lambda(i,j),n(i,j),B(i,j)] = calculate_properties_fit_simple(y(i,1),y(i,2),y(i,3),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin,Rox,Rad,L);
                [Rbulk(i,j),Rslip(i,j),vav(i,j),vslip(i,j),P(i,j), Q, Qtest(i,j), res, rexp(i,j), Pnew, A, Res2(i,j), Pressure2(i,j), Rbulk_alg(i,j) , Rslip_alg(i,j), R1(i,j), R3(i,j), Q_21, Pressure_bypass(i,j), Rslip_conQ(i,j), Rbulk_conQ(i,j),Q_21_sum, Q_byp(i,j), Q_byp_lin(i,j), Rbulk_oxy(i,j), Rslip_oxy(i,j), P2_oxy(i,j), qbulk, qslip, G(i,j)] = calculate_properties_fit_final(y(i,1),y(i,2),y(i,3),y(i,4),alpha,beta,inlet_nodes,outlet_nodes,Pin(i,1),H,D,L,z(i,:),i,j);
            end
        end
fprintf('j = %d',j)
end




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

res_final(:,1) = Rbulk;
res_final(:,2) = Rslip; 
res_final(:,3) = 1./((1./Rbulk)+(1./Rslip)); 
res_final(:,4) = Rbulk_oxy;
res_final(:,5) = Rslip_oxy; 
res_final(:,6) = 1./((1./Rbulk_oxy)+(1./Rslip_oxy)); 
res_final(:,7) = P2_oxy;
res_final(:,8) = Pressure_bypass;

%fname = sprintf('%s_resistance_analysis_%s.mat',date,id);
%save(fname,'res_final');

end




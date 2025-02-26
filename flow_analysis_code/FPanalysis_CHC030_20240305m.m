%% Flow Profiles
% Code to analyze binned data from bypass styled devices. 

%% Code and Author Details
% Data can be multiple bins decided by the user and tracked via 
% "captureBloodFlow.m". Multiple bins can be selected. 
% Manual entry of time points is required. 

% Hannah M. Szafraniec
% Living Devices Lab

%% Clear old data and initialize

close all
clear variables
clc

%% Data Import
% Insert file names WITHOUT file extension
expName = 'flowProfiles_testSP3';
oxyName = 'neofox';
sample = 'CHC030';

fprintf('Flow Profile Analysis - %s\n',expName);
data = load([expName,'.mat']);
[Time,Oxygen,TauPhaseMethod,SensorTemperature,AirPressure] = ImportOxyData([oxyName,'.csv']);

% Determine range of values (checkRange=1 if you need the first plot)
checkRange = 1;
%% Custom color codes for plotting

blue = [0 0.4470 .7410];        
orange = [.85 .325 .0980];      
yellow = [0.9290 .6940 0.1250]; 
purple = [.4940 .1840 .5560];   
green = [.4660 .6740 .1880];    
red = [228,26,28]/255;
dark_red = [173,0,0]/255;
black = [0 0 0];
g = [0.3 0.3 0.3];

%% Determining the oxygen zones
% all data
startVector = 1;
endVector = size(data.totalTime,2);

% Number of oxygen steps (up to 7 based on experimental protocol)
Psteps = 8;
%Pressure_output = 2.0;
press_step = [0,2,4,6];
O2steps = [0,2,4,6];
% Oxygen Steps
presPart = zeros(Psteps,2);

presPart(1,1) = 1; presPart(1,2) = 334; %296; %21
presPart(2,1) = 1; presPart(2,2) = 555; %364; %0
presPart(3,1) = 1; presPart(3,2) = 79; %411; %21
presPart(4,1) = 1; presPart(4,2) = 138; %468;  %2
presPart(5,1) = 1; presPart(5,2) = 163; %516; %21
presPart(6,1) = 1; presPart(6,2) = 222; %590; %4
presPart(7,1) = 1; presPart(7,2) = 240; %641; %21
presPart(8,1) = 1; presPart(8,2) = 296; %705; %6
% presPart(9,1) = 1; presPart(9,2) = 457; %750; %21
% presPart(10,1) = 1; presPart(10,2) = 523; %791; %5
% presPart(11,1) = 1; presPart(11,2) = 553; %807; %21
% presPart(12,1) = 1; presPart(12,2) = 630; %834; %7
% presPart(13,1) = 1; presPart(13,2) = 645; %850; %21
% presPart(14,1) = 1; presPart(14,2) = 663; %884; %12

%% Converting the data range

% Reading oxygen data and interpolating through missing values
formatIn = 'HH:MM:SS';
Noxy = length(Oxygen);
oxytime = datetime(Time,'InputFormat','yyyy-MM-dd HH:mm:ss.S');
flow_time = convertTo(data.frameTimeStamp_total,'datenum');
flow_time = data.frameTimeStamp_total;
day2sec = 86400;  % Seconds per day

%setup variables
pressure = data.totalPressure;
intensity = mean(data.totalIntensity,1);
flow = data.totalFlow;
intensityBin = data.intensityBin;
vel_bins = nanmean(data.velocity_bin);
velocity_bin = data.velocity_bin;
vel_bins = squeeze(vel_bins);
velocity = nanmean(vel_bins)';
Iblank = data.I_blank;
Iblood = data.I_blood;
Adj_hct = 1;%0.777251432;
hct_oxy = 0.0098*Adj_hct;
hct_deoxy = 0.0089*Adj_hct;
OD2 = data.I_blank-data.I_blood;
%OD = (log10(data.I_blank./data.I_blood));
OD = (log10(data.I_blank./data.I_blood)+log10(0.9))./hct_oxy;
OD_y = -log10(data.I_blank./(mean(data.I_blank)));
hct_oxychunk = log10(data.intensityBlank./data.intensityChunk)./hct_oxy;
hct_deoxychunk = (log10(data.intensityBlank./data.intensityChunk)+log10(0.9))./hct_deoxy;
trans = data.intensityChunk./data.intensityBlank;
%params vs. HCT

velocity_bin = data.velocity_bin;
vel_deoxy = velocity_bin(:,:,presPart(2,2)-15:presPart(2,2));
hct_deoxy = hct_deoxychunk(presPart(2,2)-15:presPart(2,2),:);

vel_oxy = velocity_bin(:,:,presPart(1,2)-10:presPart(1,2));
hct_oxy = hct_oxychunk(presPart(1,2)-10:presPart(1,2),:);

for i = 1:14
    a=(3*i-2);
    b=(3*i);
    deoxy_data(1,a:b) = hct_deoxy(i,:);
    deoxy_data(2:31,a:b) = squeeze(vel_deoxy(:,:,i))';
end

for i = 1:10
    a=(3*i-2);
    b=(3*i);
    oxy_data(1,a:b) = hct_oxy(i,:);
    oxy_data(2:31,a:b) = squeeze(vel_oxy(:,:,i))';
end


%Take last 1/4 of pressure data
deltabin = 10;
presPart(:,3) = deltabin;
presPart(:,4) = presPart(:,2) - presPart(:,3); %start vector presPart(:,4), end vector presPart(:,2)


%Binning pressure, flow, velocity data...
for i = 1:Psteps
    presBin(:,i) =  pressure(presPart(i,4):presPart(i,2),1);  %psi
end

for i = 1:Psteps
    flowBin(:,i) =  flow(presPart(i,4):presPart(i,2),1).*1000; %nL/min
end
b=1; %1 for single 20X20 um channel 2 for 40 by 20 um 
for i = 1:Psteps
    velFBin(:,i) =  (velocity(presPart(i,4):presPart(i,2),1).*400.*b.*60)./1000000;  %nL/min
end

for i = 1:Psteps
    vel(:,i) =  velocity(presPart(i,4):presPart(i,2),1);
end

vel_Bin = zeros(Psteps,30);

for i = 1:Psteps
    vel_Bin(i,:) =  nanmean(vel_bins(:,presPart(i,4):presPart(i,2)),2);
    vel_Bin2(:,:,i) = vel_bins(:,presPart(i,4):presPart(i,2));
end

vel_Bin = vel_Bin';

OD = OD';
for i = 1:Psteps
    HCT(i,:) =  OD(presPart(i,4):presPart(i,2),1);
end

HCT_avg = nanmean(HCT,1);
for i = 1:Psteps/2
HCT_oxy(i) = HCT(2*i-1);
HCT_deoxy(i) = HCT(2*i);
end

HCT_oxy = HCT_oxy';
HCT_deoxy = HCT_deoxy';

%% Only run to determine the right range of values
% Plots the complete data set to allow manual determination of individual
% sections
    
if checkRange == 1
    checkplot = figure(2);
    set(checkplot,'position',[600 100 800 800]);

    ax1 = subplot(4,1,1);
    plot(OD,'.','color',orange)
    y1 = ylabel('OD');
    %set(gca,'xtick',[])
    
    ax2 = subplot(4,1,2);
    plot(pressure,'.','color',purple);
    y2 = ylabel('Pressure (PSI)');
    set(gca,'xtick',[])
    
    ax3 = subplot(4,1,3);
    hold on
    plot(flow*1000,'.','color',green)
    plot((velocity*800*60/1000000),'.','color',blue)
    legend('Fluigent','Computer Vision')
    hold off
    y3 = ylabel('Flow Rate (nL/min)');

    ax4 = subplot(4,1,[4]);
    yyaxis right
    plot(oxytime,Oxygen,'linewidth',2)
    set( get(subplot(4,1,[4]),'YLabel'), 'String', 'Oxygen %' );
    yyaxis left
    plot(flow_time,velocity,'.','color',blue);
    xlabel('Time')
    y4 = ylabel('Velocity (um/s)');
%     ax2 = 1:1:827;
%     plot(ax2,velocity,'.','color',blue)
    ax2.XAxisLocation = 'top';
    
    ax1.FontSize = 11;
    ax2.FontSize = 11;
    ax3.FontSize = 11;
    ax4.FontSize = 11; 
    
    
    figure(3)
    yyaxis right
    plot(oxytime,Oxygen,'linewidth',3)
    yyaxis left
    plot(flow_time,velocity,'.','color',blue,'MarkerSize',10);
    xlabel('Time')
    ylabel('Velocity (um/s)');
    
    figure(4)
    plot(velocity,'.','color',blue,'MarkerSize',10);
    xlabel('Time')
    ylabel('Velocity (um/s)');
    
    
end



%% Fit Parameters
% if checkRange == 2
% fType_other = fittype('Vmax*(1-(1-wall)*abs(binNum-x)^B)',...
%                     'independent',{'binNum'},...
%                     'dependent',{'profiles'},...
%                     'coefficients',{'Vmax','wall','B','x'});
% % 
%         for i = 1:Psteps
% 
%         channel = linspace(-1,1,30);
%         binNum = 1:30;
%         channel2use = channel(binNum);
% 
%         [fitobj,gof,output] = fit(channel2use',vel_Bin(binNum,i),...
%         fType_other,'StartPoint',[vel_Bin(15,1),0.5,2,0]);
% 
%         fitted_profile(i,1:30) = fitobj(channel);
%         ci = confint(fitobj);
% % 
%         oxy_val = [21,0,21,2,21,3,21,4,21,5,21,7,21,12];
%        
%         curve_parmsA(i,1) = fitobj.wall;
%         curve_parmsA(i,2) = fitobj.x;
%         curve_parmsA(i,3) = fitobj.B;
%         curve_parmsA(i,4) = fitobj.Vmax;
%         curve_parmsA(i,5) = oxy_val(i);
%         curve_parmsA(i,6) = mean(presBin(:,i));
%         curve_parmsA(i,7) = mean(vel(:,i));
%         curve_parmsA(i,8) = mean(velFBin(:,i));
%         curve_parmsA(i,9) = mean(flowBin(:,i));
%        % curve_parms(i, 6) = vel_Bin(1,i); %top wall velocity
%        % curve_parms(i, 7) = vel_Bin(30,i); %bottom wall velocity
%  %       curve_parms(i, 7) = (vel_Bin(1,i) + vel_Bin(30,i))/2;
% %        curve_parms(i, 8) = fitobj.Vmax * fitobj.wall;
% %         curve_parms(i,9) = curve_parms(i,7)/mean(vel_Bin(15:17,i));
%          curve_parmsA(i,10) = mean(vel_Bin(:,i));
%          curve_parmsA(i,11) = max(vel_Bin(:,i));
%          curve_parmsA(i, 12) = (vel_Bin(1,i) + vel_Bin(30,i))/2;
%          curve_parmsA(i,13) = curve_parmsA(i,12)/curve_parmsA(i,11);
% %         curve_parms(i,12) = curve_parms(i,10)-curve_parms(i,11);
% %         curve_parms(i,13) = mean(presBin(:,i));
% %         curve_parms(i,14) = curve_parms(i,11)/curve_parms(i,10);
% % 
%         end
%      curve_parms(:,1) = curve_parmsA(:,13);  
%      curve_parms(:,2) = curve_parmsA(:,2); 
%      curve_parms(:,3) = curve_parmsA(:,10); 
%      curve_parms(:,4) = curve_parmsA(:,11); 
%      %curve_parms = curve_parmsA(:,1:4);
%      %curve_parms = curve_parmsA(:,1:4);
%      for i = 1:size(O2steps,2)
%          output_params(1,(4*i-3):(i*4)) = curve_parms((2*i-1),:);
%          output_params(2,(4*i-3):(i*4)) = curve_parms((2*i),:);
%      end
%      
%      
%     
% figure(3)
% hold on
% plot(channel, vel_Bin(:,1),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(1,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,2),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(2,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,4),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(4,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,6),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(6,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,8),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(8,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,10),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(10,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,12),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(12,:),'LineWidth',2,'color',orange)
% plot(channel, vel_Bin(:,14),'o','color',purple,'LineWidth',2)
% plot(channel, fitted_profile(14,:),'LineWidth',2,'color',orange)
% hold off
% xlabel('Width')
% ylabel('Velocity (um/s)');
% 
% 
% end

%% Fit Parameters
checkRange = 3;
if checkRange == 3
    fType_other = fittype('Vmax*(1-(1-wall)*abs(binNum-x)^B)',...
        'independent',{'binNum'},...
        'dependent',{'profiles'},...
        'coefficients',{'Vmax','wall','B','x'});
    %


 for j = 1:11
        vel_Bin3 = squeeze(vel_Bin2(:,j,:));

        for i = 1:Psteps

            channel = linspace(-1,1,30);
            binNum = 1:30;
            channel2use = channel(binNum);

            [fitobj,gof,output] = fit(channel2use',vel_Bin3(binNum,i),...
                fType_other,'StartPoint',[vel_Bin3(15,i),0.5,2,0]);

            fitted_profile(i,1:30) = fitobj(channel);
            ci = confint(fitobj);

            oxy_val = [21,0,21,2,21,3,21,4,21,5,21,7,21,12];

            curve_parmsA(i,1) = fitobj.wall;
            curve_parmsA(i,2) = fitobj.x;
            curve_parmsA(i,3) = fitobj.B;
            curve_parmsA(i,4) = fitobj.Vmax;
            curve_parmsA(i,5) = oxy_val(i);
            curve_parmsA(i,6) = mean(presBin(:,i));
            curve_parmsA(i,7) = mean(vel(:,i));
            curve_parmsA(i,8) = mean(velFBin(:,i));
            curve_parmsA(i,9) = mean(flowBin(:,i));
            % curve_parms(i, 6) = vel_Bin(1,i); %top wall velocity
            % curve_parms(i, 7) = vel_Bin(30,i); %bottom wall velocity
            % curve_parms(i, 7) = (vel_Bin(1,i) + vel_Bin(30,i))/2;
            % curve_parms(i, 8) = fitobj.Vmax * fitobj.wall;
            % curve_parms(i,9) = curve_parms(i,7)/mean(vel_Bin(15:17,i));
            curve_parmsA(i,10) = mean(vel_Bin3(:,i));
            curve_parmsA(i,11) = max(vel_Bin3(:,i));
            curve_parmsA(i, 12) = (vel_Bin3(1,i) + vel_Bin3(30,i))/2;
            curve_parmsA(i,13) = curve_parmsA(i,12)/curve_parmsA(i,11);
            % curve_parms(i,12) = curve_parms(i,10)-curve_parms(i,11);
            % curve_parms(i,13) = mean(presBin(:,i));
            % curve_parms(i,14) = curve_parms(i,11)/curve_parms(i,10);
            %
        end
        
        curve_parms(:,1) = curve_parmsA(:,13);
        curve_parms(:,2) = curve_parmsA(:,2);
        curve_parms(:,3) = curve_parmsA(:,10);
        curve_parms(:,4) = curve_parmsA(:,11);

        for i = 1:size(O2steps,2)
            output_params(1,(4*i-3):(i*4),j) = curve_parms((2*i-1),:);
            output_params(2,(4*i-3):(i*4),j) = curve_parms((2*i),:);
        end

   
 end
pressure = mean(presBin,1);
save( sprintf('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/input_%s_new.mat',sample), 'output_params');
save( sprintf('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/pressure_%s_new.mat',sample), 'pressure');
%csvwrite('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/input_CHC058_new.mat',output_params); 
%csvwrite('/Users/hannahszafraniec/Documents/MATLAB/network_code/input_data/pressure_CHC058_new.csv',pressure); 
end

%clearvars
%cd 'C:\Users\localuser\Documents\Hannah\jamming'
%% 
clear all
LoadPhantomLibraries();
tag = 1;
FileName = sprintf('fig5_deoxy_cut.cine');
fname =  FileName;


% %% Select roi inside channel
% [vidKLT(:,:,1), origIm] = ReadCineFileImage(fname,i,0);
%  imshow(vidKLT)
%  h = drawrectangle;
%  
% for i = 1:8000
% a = 1; %1000*(j-1);
% %b = 1000*j;
% %for i = a:(b-1)                   
%    % [vidKLT(:,:,i-a+1), origIm] = ReadCineFileImage(fname,i,0);
%     [vidKLT(:,:,a), origIm] = ReadCineFileImage(fname,i,0);
%     fprintf('i = %d\n',i);
%     vidKLT_crop = imcrop(vidKLT,h.Position);
%     %vid = im2uint8(vidKLT);
%     intensity(i,1) = mean(mean(vidKLT_crop));
% % end
% % filename = sprintf('vidKLT_%d.mat',j);
% % save(filename,'vidKLT');
% 
% end
 
% %%Select roi outside channel

% [vidKLT(:,:,1), origIm] = ReadCineFileImage(fname,1,0);
%  imshow(vidKLT)
%  h = drawrectangle;
%  
% for i = 1:8000
% a = 1; %1000*(j-1);
% %b = 1000*j;
% %for i = a:(b-1)                   
%    % [vidKLT(:,:,i-a+1), origIm] = ReadCineFileImage(fname,i,0);
%     [vidKLT(:,:,a), origIm] = ReadCineFileImage(fname,i,0);
%     fprintf('i = %d\n',i);
%     vidKLT_crop = imcrop(vidKLT,h.Position);
%     %vid = im2uint8(vidKLT);
%     intensity_blank(i,1) = mean(mean(vidKLT_crop));
% % end
% % filename = sprintf('vidKLT_%d.mat',j);
% % save(filename,'vidKLT');
% 
% end

% M = movmean(intensity,1);
% M2 = intensity_blank;
% M3 = M./(2^16);
% M3 = M3.*256;
% 
% M2 = M2./(2^16);
% M2 = M2.*256;
% %vidKLT_crop = imcrop(vidKLT,h.Position);
% 
% hct_oxy = 0.012;
% hct_deoxy = 0.0103;
% b = 0.00;
% trans = M3./M2;
% OD = -log10(trans)+b;
% HCT = (OD)./hct_deoxy;
% 
% meanHCT = mean(HCT);
% HCT_movmean = movmean(HCT,20);
% 
% save('20250123_analysis_fig5_deoxy_1.mat')


%% Get flow profile data from chunked video

clear all

LoadPhantomLibraries();
tag = 1;
FileName = sprintf('deoxy_40x.cine');
fname =  FileName;

for i = 1:4
[vidKLT(:,:,i), origIm] = ReadCineFileImage(fname,i,0);
end
%implay(vidKLT)
%%
nFrames = 4;
nFiles = 10;
conversion = 5.6/40;
nBins = 30; %set value for flow profile binning
velocity_mean = zeros(nFrames-1,nFiles);
velocity_median = zeros(nFrames-1,nFiles);
totalIntensity = zeros(nFrames-1,nFiles);

for i = 1:200
    nBins = 30;
    % read filenames
     %size(vidKLT,3);
    totalFPS(i) = 400;
    clear vidKLT
    for k = 1:4
    n = k+(nFrames*10*(i-1));
    [vidKLT(:,:,k), origIm] = ReadCineFileImage(fname,n,0);
    end
    vidKLT = im2uint8(vidKLT);



    % ROI regioning with wall removal
    if mean2(vidKLT(:,:,1)) <= 20
        objectRegion = [1 1 size(vidKLT(:,:,1),2), size(vidKLT(:,:,1),1)];
    else
        [y1,y2] = extractBoundaryMAT(vidKLT(:,:,1));
        %w = round(((y2 - y1) - (130))/2);
        width_og = y2-y1;
        blankRegion = [1 400 size(vidKLT(:,:,1),2) y1];
        %hctRegion = [z1 y1+10 z2 y2-y1-10];
        y1 = y1+0; y2 = y2-0;
        z2 = size(vidKLT(:,:,1),2);
        z1 = 1;
        %z2 = 928-z1;
        objectRegion = [z1 y1 z2 y2-y1];
        
    end

   % objectRegion = [1 175 1280 200];
    initialFrame = imcrop(vidKLT(:,:,1),objectRegion);
    
    %Signal processing code
    sig(:,1) = mean(initialFrame,1);
    sig(:,2) = movmean(sig(:,1),100,1);
    sig2 = (diff(sig(:,2)))*100;
    sig2 = movmean(sig2(:,1),100);
    [A(1),B(1)] = max(sig2);
    [A(2),B(2)] = min(sig2);
    C(1) = max(sig(:,2));
    C(2) = min(sig(:,2));
    if B(1)< B(2)
    D(1) = mean(sig(1:B(1),2));
    D(2) = mean(sig(B(1):B(2),2));
    D(3) = mean(sig(B(2):end,2));
    E(1)=B(1);E(2)=B(2);
    else
    D(1) = mean(sig(1:B(2),2));
    D(2) = mean(sig(B(2):B(1),2));
    D(3) = mean(sig(B(1):end,2));
    E(1)=B(2);E(2)=B(1);
    end
    if (z2-E(2)) < 20
       E(2) = z2;
    end
    
    if (E(1)-z1) < 20
       E(1) = z1;
    end

    n = 1;
    
   %initialFrame = imcrop(vidKLT(:,:,1),objectRegion);
    % Compute channel intensity and blank intensity for HCT estimates
    blank = imcrop(vidKLT(:,:,1),blankRegion);
    intensityBlood(i,n) = mean2(initialFrame);
    intensityBlank(i,n) = mean2(blank);

    
    % Defining Bin Widths
    width = size(initialFrame,1);
    binWidth = width/nBins;
    
    
    % Detect interest points in the object region.
    points = detectMinEigenFeatures(initialFrame, 'FilterSize',3,'MinQuality', 0.01);
    
    
    tempvel = zeros(nFrames-1,2);
    tempint = zeros(nFrames-1,1);
    tempBin = zeros(nFrames-1,nBins);

    Points_final = zeros(size(points,1),2);

    if size(points,1) > 30
        % Create a tracker object.
        tracker = vision.PointTracker('MaxBidirectionalError',3);

        % Initialize the tracker.
        initialize(tracker, points.Location, initialFrame);

        % Loop through frames and acquring velocity data
        frame = initialFrame;
        [points,validity] = step(tracker,frame);
        Points_final = points(validity,:);
        Points_final = points(:,:).*validity;
        oldValidity = validity;
        oldPoints = points;



       
         out_new = insertMarker(initialFrame,points(validity,:),'+','size',7,'Color','red');
%         out_new = initialFrame;       
%         out_new = out_new(:,:,:);
        
        while n<nFrames
            frame = imcrop(vidKLT(:,:,n+1),objectRegion);
            blank = imcrop(vidKLT(:,:,1),blankRegion);
            intensityBlood(i,n) = mean2(frame);
            intensityBlank(i,n) = mean2(blank);

            [points,validity] = step(tracker,frame);

            % Calculating the displacements
            pointDisplacement = points(:,1)-oldPoints(:,1);

            oldPoints = points;
            oldValidity = validity;
            oldValidity1 = validity;
            tempint(n) = mean2(frame);
            
            
             out = insertMarker(frame, points(validity,:),'+','size',7,'Color','blue');
%             out = frame;
             out_new = cat(4, out_new(:,:,:,:) , out(:,:,:));
            
%             %Visualize raw velocity profiles
%             vel_test(:,1) = abs(pointDisplacement)*(videoFile.fps)*conversion;
%             vel_test(:,2) = points(:,2);
%             vel_test(:,3) = points(:,1);
%             vel_out = vel_test(all(vel_test,2),:);
%             vel_out(:,2) = vel_out(:,2)*conversion;
%            
%             xvel_out = vel_out(:,3);
%             yvel_out = vel_out(:,2);
%             vvel_out = vel_out(:,1);
% 
%             figure(4+n)
%             hold on
%             plot(vel_out(:,2),vel_out(:,1),'x','Color','b');
%             hold off
%             ylabel('Vel_x (um/s)');
%             xlabel('Width (u)');
%             ylim([300 700])
%             xlim([0 21])



            %Track binned velocities
            for b=1:nBins
                tempBin2(n,b) = mean(nonzeros(pointDisplacement(validity &...
                oldValidity & (points(:,2) >= (b-1)*binWidth) &...
                (points(:,2) <= b*binWidth))));
            end
            
            for b=1:nBins
                objectRegion2 = [1 (b-1)*binWidth size(vidKLT(:,:,1),2) b*binWidth-(b-1)*binWidth];
                frame_new = imcrop(frame,objectRegion2);
                intensityBin2(b,n) = mean2(frame_new);
            end

            Points_final2 = points(:,:).*oldValidity1;

            %Signal processing code
            sig(:,1) = mean(initialFrame,1);
            sig(:,2) = movmean(sig(:,1),100);
            sig2 = (diff(sig(:,2)))*100;
            sig2 = movmean(sig2(:,1),100);
            [A(1),B(1)] = max(sig2);
            [A(2),B(2)] = min(sig2);
            C(1) = max(sig(:,2));
            C(2) = min(sig(:,2));
            if B(1)< B(2)
            D(1) = mean(sig(1:B(1),2));
            D(2) = mean(sig(B(1):B(2),2));
            D(3) = mean(sig(B(2):end,2));
            E(1)=B(1);E(2)=B(2);
            else
            D(1) = mean(sig(1:B(2),2));
            D(2) = mean(sig(B(2):B(1),2));
            D(3) = mean(sig(B(1):end,2));
            E(1)=B(2);E(2)=B(1);
            end
            
            if (abs(D(1)-D(2))) < 1
                E(1) = z1;
            end

            if (abs(D(3)-D(2))) < 1
                E(2) = z2;
            end

%             if (z2-E(2)) < 30
%                E(2) = z2;
%             end
%             
%             if (E(1)-z1) < 30
%                E(1) = z1;
%             end
            
            D(1) = mean(sig(1:E(1),2));
            D(2) = mean(sig(E(1):E(2),2));
            D(3) = mean(sig(E(2):end,2));        
            
            if E(1) == z1 
                D(1) = 0;
            end
            
            if E(2) == z2 
                D(3) = 0;
            end

            
            

            
%              E(1) = 213;   
%              E(2) = 600;
            
            ix(1) = 1;
            if mean2(vidKLT(:,:,1)) <= 20
                ix(2) = 50;
                ix(3) = 100;
                ix(4) = size(Points_final,1);
            else
            ix(2) = find(Points_final(:,1) > E(1), 1, 'first');
            ix(3) = find(Points_final(:,1) < E(2), 1, 'last');
            ix(4) = size(Points_final,1);
            end


            for k = 1:3
            
            pointDisplacement2 = pointDisplacement(ix(k):ix(k+1),1);
            validity2 = validity(ix(k):ix(k+1),1);
            oldValidity2 = oldValidity(ix(k):ix(k+1),1);
            points2 = points(ix(k):ix(k+1),:);
        
            for b=1:nBins
        %     tempBin3(k,b) = mean(nonzeros(pointDisplacement2(validity2 &...
        %     oldValidity2 & (points2(:,2) >= (b-1)*binWidth) &...
        %     (points2(:,2) <= b*binWidth))));   
            tempBin3(k,b,n) = mean(nonzeros(pointDisplacement2(validity2 &...
            oldValidity2 & (points2(:,2) >= (b-1)*binWidth) &...
            (points2(:,2) <= b*binWidth))));
            end

            tempBin3(isnan(tempBin3))=0;
            input = abs(squeeze(tempBin3(k,:,n)))*totalFPS(i)*conversion;
%             fType_other = fittype('Vmax*(1-(1-wall)*abs(binNum-x)^B)',...
%                     'independent',{'binNum'},...
%                     'dependent',{'profiles'},...
%                     'coefficients',{'Vmax','wall','B','x'});
%             channel = linspace(-1,1,30);
%             binNum = 1:30;
%             channel2use = channel(binNum);
%             
%             [fitobj,gof,output] = fit(channel2use',input(binNum)',...
%             fType_other,'StartPoint',[max(input),0.5,2,0]);
% 
%             curve_parms(k,1,n) = fitobj.wall;
%             curve_parms(k,2,n) = fitobj.x;
%             curve_parms(k,3,n) = fitobj.B;
%             curve_parms(k,4,n) = fitobj.Vmax;

            end
            
            data1(:,:,n) = abs((tempBin3(:,:,n))*totalFPS(i)*conversion);
            
            % Iterating the loop
            n=n+1;
            intD(:,n) = D';
        end

           data2 = squeeze(mean(data1,3));
           intD2 = mean(intD,2);
           width_bin(1) = E(1)-z1;
           width_bin(2) = E(2)-E(1);
           width_bin(3) = z2 - E(2);
           %curve_parms2 = squeeze(mean(curve_parms,3));
           index = [i;i;i];
           xdat = horzcat(index,intD2);
          % xdat = horzcat(xdat,(width./1023)');
          % xdat = horzcat(xdat,curve_parms2);
          % xdat = horzcat(xdat,data2);

%         for k = 1:3
%            
%            data_frame 
%         end

           I_blood(i) = mean(intensityBlood(i,:),2);
           I_blank(i) = mean(intensityBlank(i,:),2);
        
    end
    
% ydat = vertcat(ydat,xdat);
    
    tempBin3 = squeeze(mean(tempBin3,3));


% Assign output variables

    velocity_bin(:,:,i) = abs(tempBin2)*totalFPS(i)*conversion;
    velmean = squeeze(mean(velocity_bin,1));
    velmean2 = mean(velmean,1);
    velocity_bin2(:,:,i) = abs(tempBin3)*totalFPS(i)*conversion;
    totalIntensity(:,i) = tempint;
    intensityBin(:,i) = mean(intensityBin2,2);
    intensityChunk(i,:) = D;
    OD = log10(I_blank./I_blood)./1;
    trans = I_blood./I_blank;
    fprintf('i = %d\n',i);
   
end



    figure(i)
    hold on
    plot(velocity_bin2(1,:,i),'x','Color','b','MarkerSize',12,"LineWidth",4);
    plot(velocity_bin2(2,:,i),'x','Color','g','MarkerSize',12,"LineWidth",4);
    plot(velocity_bin2(3,:,i),'x','Color','r','MarkerSize',12,"LineWidth",4);
    hold off
    ylabel('Vel_x (um/s)');
    xlabel('Width (u)');
    ylim([300 500])
    xlim([0 30])




 figure(i+1)
 imshow(vidKLT(:,:,1));
 hold on
 width = size(vidKLT(:,:,1),2);
 height = size(vidKLT(:,:,1),1);
 plot([0 E(1)], [y1 y1],'color','blue',"LineWidth",3);
 plot([0 E(1)], [y2 y2],'color','blue',"LineWidth",3);
 plot([E(1) E(1)], [y1 y2],'color','green',"LineWidth",3);
 plot([E(2) E(2)], [y1 y2],'color','green',"LineWidth",3);
 plot([E(1) E(2)], [y1 y1],'color','green',"LineWidth",3);
 plot([E(1) E(2)], [y2 y2],'color','green',"LineWidth",3);
 plot([E(2) width], [y1 y1],'color','blue',"LineWidth",3);
 plot([E(2) width], [y2 y2],'color','blue',"LineWidth",3);
 hold off
 
 
 

%% analyze profiles and HCT
deoxyHCT = 0.0103;
hct_deoxychunk = -log10(intensityChunk./intensityBlank)/deoxyHCT;

vel_deoxy = velocity_bin2(:,:,:);
hct_deoxy = hct_deoxychunk(:,:);

for i = 1:200
    a=(3*i-2);
    b=(3*i);
    deoxy_data(1,a:b) = hct_deoxy(i,:);
    deoxy_data(2:31,a:b) = squeeze(vel_deoxy(:,:,i))';
end

deoxy_data(32,:) = max(deoxy_data);
deoxy_data(33,:) = (deoxy_data(2,:)+deoxy_data(31,:))/2;
deoxy_data(34,:) = deoxy_data(33,:)./deoxy_data(32,:);
meanv = mean(deoxy_data(2:31,:));

data_out(:,1) = deoxy_data(1,:)';
data_out(:,2) = deoxy_data(34,:)';
data_out(:,3) = deoxy_data(32,:)';
data_out(:,4) = deoxy_data(33,:)';
data_out(:,5) = meanv';
data_out(:,6) = data_out(:,5)-data_out(:,4);






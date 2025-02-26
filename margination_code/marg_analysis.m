%%Margination analysis 
%Created 10/26/2023 H. Szafraniec
%Code used to measure flourescent signal from stained RBCs
% flowing in microfluidic channel

%import .czi image stack?
clear all
tiff_info = imfinfo('oxy_flow_10x_rf exp2.czi - C=0.tif'); % return tiff structure, one element per image
tiff_stack = imread('oxy_flow_10x_rf exp2.czi - C=0.tif', 1)*1; % read in first image
%J = imrotate(tiff_stack,4.3);
%I = imcrop(J,[20 246 450 110]);
%tiff_stack = I;
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('oxy_flow_10x_rf exp2.czi - C=0.tif', ii)*1;
%     J = imrotate(temp_tiff,4.3);
%     I = imcrop(J,[20 246 450 110]);
    tiff_stack = cat(3 , temp_tiff, tiff_stack);
end

% %bfr = BioformatsImage('oxy_flow_10x_rf exp2.czi');
% %stack = zeros(512,512,200);
% % for i = 1:200
% %     I = getPlane(bfr, 1, 1, i)*200;
% %     I = imadjust(I);
% % 
% %     I = im2uint8(I);
% %     %I = imadjust(I);
% %     minval(i,1) = min(I(:));
% %     maxval(i,1) = max(I(:));
% % end
%   %a = double(max(minval))/double(max(maxval));
%   %b = 1;
%   %I = imadjust(I,[0 1]);
% %   imshow(I)
%   %histogram(I,30)
% 
% % J = imrotate(I,3.5);
% % J = imcrop(J,[40 40 450 450]);
% % figure(1)
% % imshow(J,[]);
% % figure(2)
% % imshow(I,[]);
% 
% 
%  for i = 1:200
%   I = getPlane(bfr, 1, 1, i);
%   J = imrotate(I,3.5);
%   J = imcrop(J,[40 40 450 450])*200;
%   J = imcrop(J,[0 350 450 70]);
% 
%   
%   B = im2uint8(J);
%   A = im2uint8(J);
%   C = imadjust(A);
%   mask(:,:,i) = B;
%   stack(:,:,i) = A;
%   stack3(:,:,i) = C;
% 
%  end
%  
%  
%  
% 
%find mask
%Frame 59 is background

%stack2 = (stack - mask(:,:,59));
for i = 1:size(tiff_info, 1)
sum1(:,i) = sum(tiff_stack(:,:,i),2);
end


sum2 = mean(sum1,2);
sum2 = sum2./(max(sum2));
sum2 = movmean(sum2,5);
%stack2 = imadjust(stack2);
%implay(stack3);
%straighten image
figure(1)
hold on
% for i = 1:200
% plot(sum1(:,i))
% end
plot(sum2(:,:),"--","LineWidth",7)
hold off

x = linspace(-1,1,size(sum2,1))';


%%deoxy

%import .czi image stack?
clear all
tiff_info = imfinfo('deoxy_flow_10x_rf exp2.czi - C=0.tif'); % return tiff structure, one element per image
tiff_stack = imread('deoxy_flow_10x_rf exp2.czi - C=0.tif', 1)*1; % read in first image
%J = imrotate(tiff_stack,4.3);
%I = imcrop(J,[20 246 450 110]);
%tiff_stack = I;
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('deoxy_flow_10x_rf exp2.czi - C=0.tif', ii)*1;
%     J = imrotate(temp_tiff,4.3);
%     I = imcrop(J,[20 246 450 110]);
    tiff_stack = cat(3 , temp_tiff, tiff_stack);
end

% %bfr = BioformatsImage('oxy_flow_10x_rf exp2.czi');
% %stack = zeros(512,512,200);
% % for i = 1:200
% %     I = getPlane(bfr, 1, 1, i)*200;
% %     I = imadjust(I);
% % 
% %     I = im2uint8(I);
% %     %I = imadjust(I);
% %     minval(i,1) = min(I(:));
% %     maxval(i,1) = max(I(:));
% % end
%   %a = double(max(minval))/double(max(maxval));
%   %b = 1;
%   %I = imadjust(I,[0 1]);
% %   imshow(I)
%   %histogram(I,30)
% 
% % J = imrotate(I,3.5);
% % J = imcrop(J,[40 40 450 450]);
% % figure(1)
% % imshow(J,[]);
% % figure(2)
% % imshow(I,[]);
% 
% 
%  for i = 1:200
%   I = getPlane(bfr, 1, 1, i);
%   J = imrotate(I,3.5);
%   J = imcrop(J,[40 40 450 450])*200;
%   J = imcrop(J,[0 350 450 70]);
% 
%   
%   B = im2uint8(J);
%   A = im2uint8(J);
%   C = imadjust(A);
%   mask(:,:,i) = B;
%   stack(:,:,i) = A;
%   stack3(:,:,i) = C;
% 
%  end
%  
%  
%  
% 
%find mask
%Frame 59 is background

%stack2 = (stack - mask(:,:,59));
for i = 1:size(tiff_info, 1)
sum1(:,i) = sum(tiff_stack(:,:,i),2);
end


sum2 = mean(sum1,2);
sum2 = sum2./(max(sum2));
sum2 = movmean(sum2,5);
%stack2 = imadjust(stack2);
%implay(stack3);
%straighten image
figure(1)
hold on
% for i = 1:200
% plot(sum1(:,i))
% end
plot(sum2(:,:),"--","LineWidth",7)
hold off

x = linspace(-1,1,size(sum2,1))';





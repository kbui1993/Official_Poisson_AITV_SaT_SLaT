% This script performs binary segmentation on a retina vessel image
% corrupted by Poisson noise and Gaussian blur.

%% load images
imagefiles = dir('Images/retina/*.gif');
nfiles = length(imagefiles);


%% Load image and modify
ii = 10; %pick a number from 1 to 20


cdata = double(imread(strcat('Images/retina/',imagefiles(ii).name)));
cdata(cdata==0) = 200;
m = (cdata==255);

%% corrupt image with gaussian blur and Poisson noise and scale it to [0,1]
rng(1234); %set seed

%create Gaussian blur and blur the image
g=fspecial('gaussian', [10 10], 2);
cdata_blurry = myconv(cdata, g);

cdata_noisy = poissrnd(double(cdata_blurry/2)); %divide peak value by 2 and then add noise
cdata_noisy = cdata_noisy/max(cdata_noisy(:)); %scale to [0,1]


%% run Poisson AITV SaT for binary segmentation
[~, idx] = Deblur_Poisson_L1mL2_2Stage(cdata_noisy, g, 22.5, 0.25, 0.8, 1, 2);
[l1ml2_dice, amax] = max([dice(double(idx==1), double(m)), dice(double(idx==2), double(m))], [], 'linear');

figure; subplot(1,3,1); imagesc(cdata); axis off; axis image; colormap gray; title('Original');
subplot(1,3,2); imagesc(cdata_noisy); axis off; axis image; colormap gray; title('Noisy');
subplot(1,3,3); imagesc(idx == amax); axis off; axis image; colormap gray; title('Segmentation');
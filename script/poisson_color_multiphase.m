% This script performs multiphase segmentation on a color image
% corrupted by Poisson noise.

%% load images
imagefiles = dir('Images/color_image/*.mat');
nfiles = length(imagefiles);

%% choose image and modify
ii = 7; %choose from 1 to 10
rng(1234);
a = struct2cell(load(strcat('Images/color_image/',imagefiles(ii).name)));
LabelMap = a{1};
max_val = max(unique(LabelMap));
rgbImage = ind2rgb(LabelMap, jet(double(max_val)));

k = poissrnd(double(uint8(rgbImage*10)));
k1 = k/10;

%% run AITV SLaT
[SLAT_AITV_result,idx1] = Poisson_L1mL2_SLaT(k1, 1.5, 0.05, 0.6, 6);
AITV_psnr = psnr(SLAT_AITV_result, rgbImage);

%% plot image
figure; 
subplot(1, 3, 1); imagesc(rgbImage); axis off; axis image; title('Original');
subplot(1, 3, 2); imagesc(k1); axis off; axis image; title('Noisy');
subplot(1, 3, 3); imagesc(SLAT_AITV_result); axis off; axis image; title('AITV SLaT');


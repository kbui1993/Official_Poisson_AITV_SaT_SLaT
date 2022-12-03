% This script performs multiphase image segmentation on a brain image
% corrupted by Poisson noise.

%% load images
imagefiles = dir('Images/retina/*.gif');
nfiles = length(imagefiles);

load('Images/brain_image/result.mat');
load('Images/brain_image/result2.mat');
load('Images/brain_image/result3.mat');
load('Images/brain_image/result4.mat');

images = {};
images{1} = image1;
images{2} = image2;
images{3} = image3;
images{4} = image4;

%% pick an image and modify
ii=2;
rng(1234);
cdata = images{ii};
cdata(cdata==0) = 10;

% get ground truth of csf, wm, and gm
csf_m = double(cdata == 48);
white_matter_m = double(cdata==154);
grey_matter_m = double(cdata==106);

% add noise and rescale
cdata_noisy = poissrnd(double(cdata/2));
cdata_noisy = cdata_noisy/max(cdata_noisy(:));

%% Segment

% perform AITV Poisson SaT
[~, idx] = Poisson_L1mL2_2Stage(cdata_noisy, 15.0, 1, 0.6, 1, 4);

% compute dice for each region
[l1l2_csf_dice, csf_idx] = max([dice(double(idx==1), double(csf_m)), dice(double(idx==2), double(csf_m)), ...
    dice(double(idx==3), double(csf_m)), dice(double(idx==4), double(csf_m))], [], 'linear');
[l1l2_wm_dice, wm_idx] = max([dice(double(idx==1), double(white_matter_m)), dice(double(idx==2), double(white_matter_m)), ...
    dice(double(idx==3), double(white_matter_m)), dice(double(idx==4), double(white_matter_m))], [], 'linear');
[l1l2_gm_dice, gm_idx] = max([dice(double(idx==1), double(grey_matter_m)), dice(double(idx==2), double(grey_matter_m)), ...
    dice(double(idx==3), double(grey_matter_m)), dice(double(idx==4), double(grey_matter_m))], [], 'linear');

% reconstruct the image
l1l2_recon_image = ones(size(cdata_noisy))*10;
l1l2_recon_image(idx == csf_idx) = 48;
l1l2_recon_image(idx == wm_idx) = 154;
l1l2_recon_image(idx == gm_idx) = 106;

% compute psnr
l1l2_psnr = psnr(double(l1l2_recon_image), double(cdata), 154);

%% plot figure
figure; subplot(1,3,1); imagesc(cdata); axis off; axis image; colormap gray; title('Original');
subplot(1,3,2); imagesc(cdata_noisy); axis off; axis image; colormap gray; title('Noisy');
subplot(1,3,3); imagesc(l1l2_recon_image); axis off; axis image; colormap gray; title('Segmentation');


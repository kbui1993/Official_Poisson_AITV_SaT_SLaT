% This script performs multiphase image segmentation on a brain image
% corrupted by Poisson noise.

imagefiles = dir('Images/brain_image/*.mat');
nfiles = length(imagefiles);

%% pick an image and modify
ii=15;
rng(1234);
a = struct2cell(load(strcat('Images/brain_image/',imagefiles(ii).name)));
cdata = a{1};
cdata(cdata==0) = 10;

% get ground truth of csf, wm, and gm
csf_m = double(cdata == 48);
white_matter_m = double(cdata==154);
grey_matter_m = double(cdata==106);

% add noise and rescale
cdata_noisy = poissrnd(double(cdata/8));
cdata_noisy = cdata_noisy/max(cdata_noisy(:));


%% Segment

% perform AITV Poisson SaT
[~, idx] = Poisson_L1mL2_2Stage(cdata_noisy, 4., 1., 0.6, 1, 4);

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

%% plot figure
figure; subplot(1,3,1); imagesc(cdata); axis off; axis image; colormap gray; title('Original');
subplot(1,3,2); imagesc(cdata_noisy); axis off; axis image; colormap gray; title('Noisy');
subplot(1,3,3); imagesc(l1l2_recon_image); axis off; axis image; colormap gray; title('Segmentation');
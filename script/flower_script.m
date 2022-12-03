% This script performs color image segmentation on an image of a flower
% with k=8 regions.

%% load image
flower = imread("flower.jpg");

%% add noise to image
rng(1234);
flower_noisy = poissrnd(double(flower));
flower_noisy = rescale_color_image(flower_noisy);

%% perform piecewise image reconstruction
AITV_SLAT_result = Poisson_L1mL2_SLaT(flower_noisy, 10, 1, 0.8, 8);
l1ml2_psnr = psnr(rescale_color_image(AITV_SLAT_result), rescale_color_image(double(flower)));

%% plot figure
figure; subplot(1,3,1); imagesc(flower); axis off; axis image; title('Original');
subplot(1,3,2); imagesc(flower_noisy); axis off; axis image; title('Noisy');
subplot(1,3,3); imagesc(AITV_SLAT_result); axis off; axis image; title('AITV SLaT')
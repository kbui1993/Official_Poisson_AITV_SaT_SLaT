% This script performs color image segmentation on an image of a tree
% with k=8 regions.
% Note: result differs slightly from the publication.

%% load image
tree = imread("tree.jpg");

%% add noise to image
rng(1234);
tree_noisy = poissrnd(double(tree));
tree_noisy = rescale_color_image(tree_noisy);

%% perform piecewise image reconstruction
AITV_SLAT_result = Poisson_L1mL2_SLaT(tree_noisy, 8, 1, 0.8, 8);
l1ml2_psnr = psnr(rescale_color_image(AITV_SLAT_result), rescale_color_image(double(tree)));

%% plot figure
figure; subplot(1,3,1); imagesc(tree); axis off; axis image; title('Original');
subplot(1,3,2); imagesc(tree_noisy); axis off; axis image; title('Noisy');
subplot(1,3,3); imagesc(AITV_SLAT_result); axis off; axis image; title('AITV SLaT');
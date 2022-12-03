% This script performs color image segmentation on an image of a shoe
% with k=8 regions.

%% load image
shoe = imread("shoe.jpg");

%% add noise to image
rng(1234);
shoe_noisy = poissrnd(double(shoe));
shoe_noisy = rescale_color_image(shoe_noisy);

%% perform piecewise image reconstruction
AITV_SLAT_result = Poisson_L1mL2_SLaT(shoe_noisy, 8.5, 1.0, 0.8, 8);
l1ml2_psnr = psnr(rescale_color_image(AITV_SLAT_result), rescale_color_image(double(shoe)));

%% plot figure
subplot(2,3,1);
imagesc(shoe); axis off; axis image; title('Original');
% This line is to select the patch to be zoomed in :
rectangle('Position',[300 175 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height

% Draw zoomed in portion
hold on

% These two lines are  for the position of zoomed part :
subplot(2,3,4); imagesc(shoe(175:225, 300:350,:));hold on; axis off; axis square;

subplot(2,3,2); imagesc(shoe_noisy); axis off; axis image; title('Noisy');
% This line is to select the patch to be zoomed in :
rectangle('Position',[300 175 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height

% Draw zoomed in portion
hold on

subplot(2,3,5); imagesc(shoe_noisy(175:225, 300:350,:));hold on; axis off; axis square;

subplot(2,3,3);
AITV_SLAT_result = rescale_color_image(AITV_SLAT_result); 
imagesc(AITV_SLAT_result); axis off; axis image;  title('AITV SLaT');

% This line is to select the patch to be zoomed in :
rectangle('Position',[300 175 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height


% % This line connect the patch and zoomed part : 
% line([150 120], [170 202], 'Color', 'red') % right lower corner of the rectagle is (103,80)


% Draw zoomed in portion
hold on
subplot(2,3,6); imagesc(AITV_SLAT_result(175:225, 300:350,:));hold on; axis off; axis square;
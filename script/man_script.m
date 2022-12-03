% This script performs color image segmentation on an image of a man
% with k=6 regions.

%% load image
man = imread("man.jpg");

%% add noise to image
rng(1234);
man_noisy = poissrnd(double(man));
man_noisy = rescale_color_image(man_noisy);

%% perform piecewise image reconstruction
AITV_SLAT_result = Poisson_L1mL2_SLaT(man_noisy, 15, 1.0, 0.8, 6);
l1ml2_psnr = psnr(rescale_color_image(AITV_SLAT_result), rescale_color_image(double(man)));

%% plot figure
figure;
subplot(2,3,1);
imagesc(man); axis off; axis image; title('Original');
% This line is to select the patch to be zoomed in :
rectangle('Position',[75 75 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height

% Draw zoomed in portion
hold on

% These two lines are  for the position of zoomed part :
subplot(2,3,4); imagesc(man(75:125, 75:125,:));hold on; axis off; axis square;

subplot(2,3,2); imagesc(man_noisy); axis off; axis image; title('Noisy');
% This line is to select the patch to be zoomed in :
rectangle('Position',[75 75 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height

% Draw zoomed in portion
hold on

subplot(2,3,5); imagesc(man_noisy(75:125, 75:125,:));hold on; axis off; axis square;

subplot(2,3,3);
AITV_SLAT_result = rescale_color_image(AITV_SLAT_result);
imagesc(AITV_SLAT_result); axis off; axis image; title('AITV SLaT');

% This line is to select the patch to be zoomed in :
rectangle('Position',[75 75 50 50], 'EdgeColor','r', 'LineWidth', 3) % left upper corner of the rectagle is (52,53), then add wide and height


% % This line connect the patch and zoomed part : 
% line([150 120], [170 202], 'Color', 'red') % right lower corner of the rectagle is (103,80)


% Draw zoomed in portion
hold on
subplot(2,3,6); imagesc(AITV_SLAT_result(75:125, 75:125,:));hold on; axis off; axis square;

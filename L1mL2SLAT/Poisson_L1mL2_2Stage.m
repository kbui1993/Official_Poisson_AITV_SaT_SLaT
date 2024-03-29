%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function performs AITV two-stage segmentation for grayscale images.
%Input:
%   f: image
%   lambda: weighing parameter for fidelity term
%   mu: weigh parameter for smoothing term
%   alpha: sparsity parameter for L1-\alpha L2 term of gradient
%   beta: penalty parameter for ADMM
%   k: number of regions/clusters
%Output:
%   result: segmented image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,idx] = Poisson_L1mL2_2Stage(f, lambda, mu, alpha, beta, k)

    %obtain size of f
    [m,n] = size(f);
    
    %% Two-Stage Segmentation
    %stage one: smooth the image
    u = Poisson_L1mL2smooth(f, lambda, mu, alpha, beta);
    
    %stage two: perform k-means
    u_vector = reshape(u, m*n,1);
    idx = kmeans(u_vector, k, 'Replicates',5);
    idx = reshape(idx, m,n);
    
    %% postprocessing step
    %compute mean of each part and construct the piecewise constant image
    mean_c = zeros(k,1);
    result = zeros(m,n);
    for i=1:k
        region = (i==idx).*u;
        mean_c(i) = sum(region(:))/(sum(region(:)>0));
        result = result+(i==idx)*mean_c(i);
    end
end
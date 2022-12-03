%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function solves the following nonconvex variant of the Mumford-Shah
%model with AITV regularization:
%
%   min \lambda \langle Au - f log Au \rangle + \frac{\mu}{2} \|\nabla u\|_F^2 +
%   \|\nabla u\|_1 - \alpha \|\nabla u\|_{2,1}
%
%Input:
%   f: image
%   A: deblurring operator
%   lambda: weighing parameter for fidelity term
%   mu: weighing parameter for smoothing term
%   alpha: sparsity parameter for L1-\alpha L2 term of gradient
%   beta: penalty parameter for ADMM
%
%Output:
%   u: solution/smoothed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u]= Deblur_Poisson_L1mL2smooth(f, A, lambda, mu, alpha, beta)
    
    %penalty parameter multiplier
    rho = 1.25;
    
    %obtain dimension of image
    [rows,cols] = size(f);
    
    %preintialize variable to store past u
    u0 = ones(rows,cols);
    
    %preinitialize u
    u= u0;

    %preinitialize v
    v=u0;
    
    %preinitialize w variables
    wx = u0;
    wy = u0;
    
    %preinitialize dual variable
    z = v;
    zx = u0;
    zy = u0;
    
    %build kernel: use the fft algorithm (5-pt stencil)
    uker = zeros(rows,cols);
    uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
    
    %refit blurring operator and shift it
    [xLen_flt, yLen_flt] = size(A);
    ope_blur=zeros(rows,cols);
    ope_blur(1:xLen_flt,1:yLen_flt)=A;
    
    xLen_flt_1=floor(xLen_flt/2);yLen_flt_1=floor(yLen_flt/2);
    ope_blur_1=padarray(ope_blur,[rows,cols],'circular','pre');
    ope_blur_1=ope_blur_1(xLen_flt_1+1:rows+xLen_flt_1,yLen_flt_1+1:cols+yLen_flt_1);
    
    %fourier transform of blurring operator
    FA = fft2(ope_blur_1);

    %compute Au
    Au = ifft2(FA.*fft2(u));
    
    for i=1:300

 %store past u
        u0 = u;
        
        %left-hand side of optimality eqn of u
        new_uker = beta*conj(FA).*FA+(mu+beta)*fft2(uker);
        
        %right-hand side of optimality eqn of u
        rhs1 = beta*v-z;
        rhs2 = beta*Dxt(wx)-Dxt(zx)+beta*Dyt(wy)-Dyt(zy);
        
        %solve u-subproblem
        u = ifft2((conj(FA).*fft2(rhs1)+fft2(rhs2))./new_uker);

        %compute Au
        Au = ifft2(FA.*fft2(u));
        
        %compute relative err
        err=norm(u-u0,'fro')/norm(u, 'fro');
        
        if mod(i,10)==0
            disp(['iterations: ' num2str(i) '!  ' 'error is:   ' num2str(err)]);
        end
        
        % check the stopping criterion
        if err<10^(-4)
            break;
        end

        %solve v-subproblem
        v = (sqrt((lambda-z-beta*Au).^2+4*beta*lambda*f)-(lambda-z-beta*Au))/(2*beta);

        
        %solve w-subproblem
        temp1 = Dx(u)+zx/beta;
        temp2 = Dy(u)+zy/beta;
        
        temp1 = reshape(temp1, rows*cols,1);
        temp2 = reshape(temp2, rows*cols,1);
        
        temp = [temp1, temp2];
        temp = shrinkL12(temp,1/beta, alpha);
        wx = temp(:,1);
        wy = temp(:,2);
        wx = reshape(wx, rows,cols);
        wy = reshape(wy, rows,cols);
        
        %update dual variables
        zx = zx+beta*(Dx(u)-wx);
        zy = zy+beta*(Dy(u)-wy);
        z = z+beta*(Au-v);
        
        %update ADMM penalty parameter
        beta = beta*rho;
    end


end

function x = shrinkL12(y,lambda,alpha)
    %this function applies the proximal operator of L1-alpha L2 to each
    %row vector
    
    %initialize solution as zero vector
    x = zeros(size(y));
    
    %obtain the indices of the max entries of each row vector
    [max_y, idx_y] = max(abs(y'));
    max_y = max_y';
    idx_y = idx_y';
    new_idx_y = sub2ind(size(y), (1:size(y,1))',idx_y);
    
    %compute new row vectors when max value of each row vector is greater
    %than lambda
    case1_idx = max_y > lambda;
    
    case1_result = max(abs(y(case1_idx,:))-lambda,0).*sign(y(case1_idx,:));
    norm_case1_result = sqrt(sum(case1_result.^2,2));
    x(case1_idx,:) =((norm_case1_result+alpha*lambda)./norm_case1_result).*case1_result;
    
    %compute one-sparse vector when max value of each row vector is less
    %than or equal to lambda and above (1-alpha)*lambda
    case2_idx = logical((max_y<=lambda).*(max_y>=(1-alpha)*lambda));
    
    x(new_idx_y(case2_idx)) = (max_y(case2_idx)+(alpha-1)*lambda).*sign(y(new_idx_y(case2_idx)));
    
end


    
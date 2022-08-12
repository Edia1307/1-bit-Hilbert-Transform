%% BIHT-l1 Subfunction
function [ x_Norm ] = sub_BIHT_l1(Phi,y,K)
    [~,N] = size(Phi);
    sgn =  @(in) sign(in+eps);
    A = @(in) sgn(Phi*in);  % complex sign function
    maxiter = 200;
    htol = 0;
    x = zeros(N,1);
    hd = Inf;
    ii=0;%
    while (htol < hd)&&(ii < maxiter)
        % Get gradient
        g = Phi'*(A(x) - y);
        % Step
        a = x - g;
        % Best K-term (threshold)
        [~, aidx] = sort(abs(a), 'descend');
        a(aidx(K+1:end)) = 0;
    % 	temp_col(:,ii+1) = sort(aidx(1:K),'ascend'); 
        % Update x
        x = a;
        % Measure hammning distance to original 1bit measurements
        hd = nnz(y - A(x));
        ii = ii+1;
    end
    % Now project to sphere
    x_Norm = x/norm(x);
end
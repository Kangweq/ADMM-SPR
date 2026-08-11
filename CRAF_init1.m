function x_est = CRAF_init1(y, A, k)
    [m, n] = size(A);
    y_abs2 = y.^2;
    
   
    lambda_neg = -3;
    lambda_pos = 1;
    
    
    r_hat = sqrt(mean(y.^2));
    I_neg = find(y.^2 <= r_hat^2 / 2);
    I_pos = find(y.^2 >= r_hat^2 / 2);
    
    
    zeta = ((y_abs2)' * (abs(A).^2))' / m;
    [~, supp_est] = maxk(zeta, k);   
    
   
    A_supp = A(:, supp_est);
    A_neg = A_supp(I_neg, :);   
    A_pos = A_supp(I_pos, :);   
    len_neg = length(I_neg);
    len_pos = length(I_pos);
    
    
    c = 5;   
    z0 = randn(k, 1);
    z0 = z0 / norm(z0);
    
    npower_iter = 100;   
    for iter = 1:npower_iter
        %  M * z0 = (λ⁻/|I⁻|) * A_neg' * (A_neg*z0) + (λ⁺/|I⁺|) * A_pos' * (A_pos*z0)
        Mz = (lambda_neg / len_neg) * (A_neg' * (A_neg * z0)) + ...
             (lambda_pos / len_pos) * (A_pos' * (A_pos * z0));
        %  (M + cI) * z0 = Mz + c*z0
        z0 = Mz + c * z0;
        z0 = z0 / norm(z0);
    end
    
    d_supp = z0;   
    
    
    d = zeros(n, 1);
    d(supp_est) = d_supp;
    x_est = r_hat * d;
end
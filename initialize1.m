function x0 = initialize1(y,A,s)
    % support recovery
    y_abs2=y.^2;
    [m,n]=size(A);
    Marg=((y_abs2)'*(abs(A).^2))'/m;%y'*A.^2;
    [~,MgS] = sort(Marg,'descend');
    S0 = MgS(1:s); %pick top s-marginals
    Shat = sort(S0); %store indices in sorted order
    A1 = A(:,Shat);
    npower_iter = 100;                          % Number of power iterations 
    phi_sq = sum(y_abs2)/m;
    phi = sqrt(phi_sq); 
    z0 = randn(s,1); z0 = z0/norm(z0,'fro');    % Initial guess 
    for tt = 1:npower_iter                     % Power iterations
        z0 = A1'*(y_abs2.* (A1*z0))/m; z0 = z0/norm(z0,'fro');
    end  
    z1=zeros(n,1);
    z1(Shat)=z0;
    x0 = phi * z1;                           % Apply scaling 
end
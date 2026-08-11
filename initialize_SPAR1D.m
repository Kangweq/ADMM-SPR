function z = initialize_SPAR1D(y,Params,Amatrix) % #ok<INUSD>

Arnorm  = sqrt(sum(abs(Amatrix).^2, 2)); % norm of rows of Amatrix
ymag    = sqrt(y);
normest = Params.m * Params.n1 / sum(sum(abs(Amatrix))) * (1 / Params.m) * sum(ymag);
ynorm   = ymag ./ (Arnorm .* normest);

%% finding largest normalized inner products
Anorm      = bsxfun(@rdivide, Amatrix, Arnorm);
[ysort, ~] = sort(ynorm, 'ascend');
ythresh    = ysort(round(Params.m / (1.2))); % 6/5 the orthogonality-promoting initialization parameter
ind        = (abs(ynorm) >= ythresh);

%% estimate the support of x
Aselect  = Anorm(ind, :);

% based on orthogonality-promoting initialization
rdata_opi= sum(abs(Aselect).^2, 1);
[~, sind_opi] = sort(rdata_opi, 'descend');
Supp_opi = sind_opi(1 : round(Params.n1 - Params.nonK));

  Asample  = Amatrix(:, Supp_opi);
Arnormx  = sqrt(sum(abs(Asample).^2, 2)); % norm of rows of Amatrix

% finding largest normalized inner products
Anormx   = bsxfun(@rdivide, Amatrix, Arnormx);
ynormx   = ymag ./ (Arnormx .* normest);
ysortx   = sort(ynormx, 'ascend');

ythreshx = ysortx(round(Params.m / (1.2))); % 6/5 the orthogonality-promoting initialization parameter
indx     = (abs(ynormx) >= ythreshx);
Aselectx = Anormx(indx, Supp_opi);
%Y=Aselectx'*Aselectx;

% [V, D] = eig(Y);  % Eigen decomposition
% [~, idx1] = max(diag(D));  % Find eigenvector corresponding to largest eigenvalue
% zk0= V(:, idx1);

zk0      = randn(Params.n1 - Params.nonK, Params.n2);
zk0      = zk0 / norm(zk0, 'fro');    % Initial guess
for t    = 1:Params.npower_iter                   % Truncated power iterations
    zk0  = Aselectx' * (Aselectx * zk0);
    zk0  = zk0 / norm(zk0, 'fro');
end

z0      = zeros(Params.n1, 1);
z0(Supp_opi) = zk0; 
%end
z       = normest * z0;
end
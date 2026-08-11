function [y_abs,A] = measure_signal(m,z,cplx_flag)
%edited 2/15/2017
n = length(z);
%% signal measurement
if cplx_flag == 0
    A  = randn(m, n); % real measurements
else
    A  = (randn(m, n) + 1i * randn(m, n)) / sqrt(2); % complex measurements
end
y = A*z; %measurements
y_abs = abs(y);
end

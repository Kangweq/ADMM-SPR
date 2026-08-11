function [z,z_ind] = generate_signal(n,K)
        z_ind = randperm(n,K); %creates sparse vector
        z = zeros(n,1);
        z(z_ind) = randn(K,1); %generate K sparse signal (n x 1)
        z = z/norm(z);   %need not be normalized     
end
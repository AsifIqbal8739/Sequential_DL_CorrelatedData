% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

%% Correlated Noise generation function for STDL (A1A2A3) 
function [Noise,Omega,Lambda,Covv] = Noise_Corr2(S,sigma)

% Defaults
    theta = 0.8;    n = S(1);   N = S(2);
    Lambda = thiss(theta,N);
    Omega = inv(AR1(n,sigma));

    Omega_sq = sqrtm(Omega);
    Lambda_sq = sqrtm(Lambda);
    Covv = sparse(kron(Lambda_sq,Omega_sq));
%     Covv = sparse(kron(Omega_sq,Lambda_sq));
    NN = randn(1,n*N);
    Noise = reshape(NN * Covv,n,N);

end

% function to make Lambda = theta^|i-j|
function tt = thiss(theta,M)
    tt = zeros(M);
    for i = 1:M
        for j = 1:M
            tt(i,j) = theta^(abs(i-j));
        end
    end
end

% Functin to make Omega = AR(1) 
function tt = AR1(m,sigma)
    a = 0.5;    
    dup = -a * diag(ones(1,m-1),1);
    ddn = -a * diag(ones(1,m-1),-1);
    dd = (1+a^2)*eye(m);  dd(1) = 1; dd(end) = 1;
    tt = 1/sigma^2 * (dd+ddn+dup);
end
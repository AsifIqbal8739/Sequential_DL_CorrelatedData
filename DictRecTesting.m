% Script to test dictionary recovery in the presence of Correlated Noise
% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

close all; clear all;
% clc; 
rng('Default');
addpath(strcat(pwd,'\ompbox10'));       % OMP box should be compiled before executing

%% Data Stuff
m = 15;     n = 40;     N = 700;       % D(m,n), Y(m,N)
K = 2;                                  % Sparsity
Sigma = 0.1;                           % Noise_Corr Setup
noIt = 11*K^2;      nTrials = 1;
[Count_K,Count_A1,Count_A2,Count_A3] = deal(zeros(nTrials,1)); 

%% Dict Creation and Signal Setup
Dict_O = normc(randn(m,n));
[Noise,Omega,Lambda,Covv] = Noise_Corr2([m,N],Sigma);

%% Trials and Signal Creation
for tr = 1:nTrials
    X = zeros(n,N);
    for i = 1:N
        y = randperm(n,K);
        X(y,i) = randn(1,K);
    end
    Y = Dict_O * X;
    Yn = Y + reshape(Noise,size(Y));
%     Yn = Y + Sigma*randn(size(Y));
    Dict = normc(Yn(:,randperm(N,n)));      % Starting Dictionary

    % Dict Recovery
    D_KSVD = K_SVD(Yn,Dict,noIt,K,1,0);
    Count_K(tr,1) = NumAtomRec(D_KSVD,Dict_O);
    D_A1 = Algo_A1(Yn,Dict,K,noIt,pinv(Omega),pinv(Lambda));
    Count_A1(tr,1) = NumAtomRec(D_A1,Dict_O);
    D_A2 = Algo_A2(Yn,Dict,K,noIt,pinv(Omega),pinv(Lambda));
    Count_A2(tr,1) = NumAtomRec(D_A2,Dict_O);
    D_A3 = Algo_A3(Yn,Dict,K,noIt,pinv(Omega),pinv(Lambda),0.2);
    Count_A3(tr,1) = NumAtomRec(D_A3,Dict_O);    
end
fprintf('Recovered Atoms by KSVD:%0.2f, A1:%0.2f, A2:%0.2f, A3:%0.2f\n',mean(Count_K),mean(Count_A1),mean(Count_A2),mean(Count_A3));

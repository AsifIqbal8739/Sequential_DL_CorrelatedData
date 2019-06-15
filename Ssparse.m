% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

%% Sparse Coding Function
function X = Ssparse(Dict,Y,Q,R,s)
N = size(Y,2);
if N <= 100
    Y_ = Y(:);
    Dict_ = normc(sparse(kron(eye(size(Y,2)),Dict)));
    LQ = sqrtm(Q);  LR = sqrtm(full(R));
    Z = sparse(kron(LQ',LR'));
    Sp = s * N;
    
    Ynew = Z * Y_;    Dnew = Z * Dict_;
    X_ = omp(Dnew,Ynew,Dnew'*Dnew,Sp);
    X = reshape(X_,size(Dict,2),N);
else % For large datasets Approximating with OMP
    X = omp(Dict'*Y,Dict'*Dict,s);
end
end
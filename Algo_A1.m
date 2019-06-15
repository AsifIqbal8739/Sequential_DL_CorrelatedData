% Code for Algorithm A1
% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

function [Dict,X] = Algo_A1(Y,Dini,s,noIt,Q,Rc)
% Y = Data Matrix
% Dini = Initial Dictionary
% s = Signal Sparsity,  noIt = Number of Dict Learning Iterations
% Q = Inverse Column Covariance, Rc = Inverse Row Covariance (Should be fullrank)

%% Parameter Setup
Comp = 0;       % Q and R matrices are to be computed from the data?
Dict = Dini;
for it = 1:noIt
    X = Ssparse(Dict,Y,Q,Rc,s);
    E = Y - Dict*X;
    if Comp == 1
        Q = pinv(cov(E'));
%         Q = inv(cov(E') + eye(size(E,1)));
    end
    LQ = sqrtm(Q);
    
    for k = 1:size(Dict,2)
        I = find(X(k,:));       
        if isempty(I)
            fprintf('Unused Atom Found!!!\n'); 
            Dict(:,k) = rep_atom(E); continue;   
        end
        Eki = E(:,I) + Dict(:,k)*X(k,I);
        if Comp == 1
            R = pinv(cov(Eki));  LR = sqrtm(R);
%             R = pinv(cov(Eki) + eye(size(Eki,1)));  LR = sqrtm(R);
        else
        R = Rc(I,I);     LR = sqrtm(R);
        end
        [U,S,V] = svds(LQ'*Eki*LR,1,'L');
        Dict(:,k) = (LQ')\U;
        X(k,I) = S*V'/LR;
        E(:,I) = Eki - Dict(:,k)*X(k,I);
    end
    Dict = normc(Dict); % Normalize Atoms
    Dinn = norm(Dict - Dini,'fro')/norm(Dini,'fro');
    fprintf('Iteration # %d, dict convergence rate : %0.4f, and E : %0.3f\n',it,Dinn,norm(Y - Dict*X,'fro')/norm(Y,'fro'));
    if Dinn <= 0.01;   break;   end
    Dict = I_clearDictionary(Dict,X,Y,1);
    Dini = Dict;
end
end

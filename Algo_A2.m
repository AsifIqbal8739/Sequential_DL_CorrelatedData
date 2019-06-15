% Code for Algorithm A2
% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

function [Dict,X] = Algo_A2(Y,Dini,s,noIt,Q,Rc)
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
    
    for k = 1:size(Dict,2)
        I = find(X(k,:));       
        if isempty(I)
            fprintf('Unused Atom Found!!!\n'); 
            Dict(:,k) = rep_atom(E); continue;   
        end
        Eki = E(:,I) + Dict(:,k)*X(k,I);
        if Comp == 1
            R = pinv(cov(Eki));
%             R = inv(cov(Eki) + eye(size(Eki,1)));  
        else
        R = Rc(I,I); 
        end
        d = Dict(:,k);
        for i =1:2
            d = Eki*R*X(k,I)';      
            d = d/sqrt(d'*Q*d);
            X(k,I) = d'*Q*Eki;
        end
        Dict(:,k) = d;      
        E(:,I) = Eki - Dict(:,k)*X(k,I);
    end
    Dict = normc(Dict); 
    Dinn = norm(Dict - Dini,'fro')/norm(Dini,'fro');
    fprintf('Iteration # %d, dict convergence rate : %0.4f, and E : %0.3f\n',it,Dinn,norm(Y - Dict*X,'fro')/norm(Y,'fro'));
    if Dinn <= 0.01;   break;   end
    Dict = I_clearDictionary(Dict,X,Y,1);
    Dini = Dict;
end
end

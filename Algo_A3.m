% Code for Algorithm A3
% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

function [Dict,X] = Algo_A3(Y,Dini,s,noIt,Q,Rc,alpha)
% Y = Data Matrix
% Dini = Initial Dictionary
% s = Signal Sparsity,  noIt = Number of Dict Learning Iterations
% Q = Inverse Column Covariance, Rc = Inverse Row Covariance (Should be fullrank)
% alpha = Smoothness Penalty Parameter

%% Parameter Setup
Comp = 0;       % Q and R matrices are to be computed from the data?
n = size(Dini,1);
t1 = eye(n)*2;    t2 = ones((n-1),1);
Omega = diag(t2,1) + diag(t2,-1) - t1; 
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
            frprintf('Unused Atom Found!!!\n'); 
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
            d = pinv((X(k,I)*R*X(k,I)')*Q + alpha*Omega)*(Q*Eki*R*X(k,I)');    %% Problem here  
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

% Written by: Asif Iqbal (aiqbal1@student.unimelb.edu.au)

%% KSVD implementation
function [D,W] = K_SVD(Y,D,noIt,s,st,d)
% st = 1st element to be optimized or not, 1 or 2 (starting atom)
if ~exist('st')  % to show the dict convergence or not.
    d = 0;  st = 1;
elseif ~exist('d')
    d = 0;
end
K = size(D,2);
D = normc(D);
for it = 1:noIt    
    Dini = D;
    W = omp(D'*Y,D'*D,s);
    %% K-SVD
    R = Y - D*W;
    for k = st:K
        I = find(W(k,:));
        if isempty(I)
            fprintf('Unused Atom Found!!!\n'); 
            Dict(:,k) = rep_atom(R); continue;   
        end
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end   
        if d == 1
            Dinn = norm(D - Dini,'fro')/norm(Dini,'fro');
            fprintf('Iteration # %d, dict convergence rate : %0.4f, and E : %0.3f\n',it,Dinn,norm(Y - D*W,'fro')/norm(Y,'fro'));
            if Dinn<=0.01;   break;   end
        end
        D = I_clearDictionary(D,W,Y,st);
end
% W = full(W);
end

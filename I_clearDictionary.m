%% Function to replace unused / similar dictionary atoms
function Dictionary = I_clearDictionary(Dictionary,W,Data,st)
if ~exist('st')  % to show the dict convergence or not.
    st = 1;
end
T2 = 0.99;
T1 = 3;
K = size(Dictionary,2);
Er = sum((Data-Dictionary*W).^2,1); % remove identical atoms
G = Dictionary'*Dictionary; G = G-diag(diag(G));
for jj = st:1:K,
    if max(G(jj,:)) > T2 || nnz(W(jj,:)) <= T1 ,
        [val,pos] = max(Er);
        Er(pos(1)) = 0;
        Dictionary(:,jj) = Data(:,pos(1))/norm(Data(:,pos(1)));
        G = Dictionary'*Dictionary; G = G-diag(diag(G));
    end
end
end
% Replace atom which isn't used
function [d] = rep_atom(Y)
    temp = sqrt(diag(Y'*Y));    % Calculating norm2 of all vectors in Y
    [a,b] = sort(temp,1,'descend');
    d = normc(Y(:,b(1)));
end
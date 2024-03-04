function adj = eignorm(mat)
% normalisation of matrix by largest eigenvalue    
    N = size(mat,2);
%     [V, D] = eig(mat);
    [U,D,V] = svd(mat);
    D2 = D./max(max(D)).*10;
    adj = (U*D2*V').*~eye(N,N);
end

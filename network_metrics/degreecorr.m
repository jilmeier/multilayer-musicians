function dc = degreecorr(W)
% compute degree correlations between several layers
% input is a W matrix with size N x N x M, where N are the no of nodes and
% M are the number of layers
M = size(W,3);

% degree per layer
degrees_layer = squeeze(sum(W,1));

% overlap degree
dc = corr(degrees_layer).*~eye(M,M);
end
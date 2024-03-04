function od = overlap_degree(W)
% compute overlap degree based on several layers
% input is a W matrix with size N x N x M, where N are the no of nodes and
% M are the number of layers

% degree per layer
degrees_layer = squeeze(sum(W,1));

% overlap degree
od = squeeze(sum(degrees_layer,2));
end
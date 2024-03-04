function H = entropylayer(W)
% compute nodal entropy based on several layers
% input is a W matrix with size N x N x M, where N are the no of nodes and
% M are the number of layers
M = size(W,3);
degrees_layer = squeeze(sum(W,1));
od = overlap_degree(W);
od_rep = repmat(od,1,M);
P = degrees_layer./od_rep; 
H = - sum(P.*log(P),2);
end
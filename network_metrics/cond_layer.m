function H = cond_layer(W)
% compute cond. probability of findings same links in several layers
% input is a W matrix with size N x N x M, where N are the no of nodes and
% M are the number of layers
M = size(W,3);
H = zeros(M,M);
for layer1 = 1:M
    for layer2 = 1:M
        if layer1 ~= layer2
        l1 = W(:,:,layer1);
        l2 = W(:,:,layer2);
        H(layer1,layer2) = sum(sum(l1.*l2))./sum(sum(l2));
        end
    end
end
end
function M = get_mutualinformation(nodeProbs, edgeProbs, numANodes, edgeEnds)

M = zeros(numANodes);

for i = 1:size(edgeEnds, 1)
    x = edgeEnds(i,1);
    y = edgeEnds(i,2);
    if x > numANodes || y > numANodes
        continue;
    end
    Px = nodeProbs(:, x);
    Py = nodeProbs(:, y);
    Pxy = edgeProbs(:,:,i);
%     M(x,y) = getEntropy(Px) + getEntropy(Py) - jointEntropy(Pxy);
    temp = log2(Pxy./(Px*Py'));
    temp(isinf(temp)) = 0;
    temp(isnan(temp)) = 0;
    M(x,y) = sum(sum(Pxy .* temp));
end

function Ke = elementStiffness(nodeCoords, D)
    gaussPoints = [-1 / sqrt(3), -1 / sqrt(3); 
                    1 / sqrt(3), -1 / sqrt(3); 
                   -1 / sqrt(3),  1 / sqrt(3); 
                    1 / sqrt(3),  1 / sqrt(3)];
    weights = [1, 1, 1, 1];
    Ke = zeros(8, 8);
    for i = 1:size(gaussPoints, 1)
        xi = gaussPoints(i, 1);
        eta = gaussPoints(i, 2);
        [B, detJ] = strainDisplacementMatrix(nodeCoords, xi, eta);
        Ke = Ke + (B' * D * B ).*( detJ * weights(i));
    end
end

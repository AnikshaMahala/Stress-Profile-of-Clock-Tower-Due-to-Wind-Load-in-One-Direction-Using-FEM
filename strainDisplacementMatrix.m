function [B, detJ] = strainDisplacementMatrix(nodeCoords, xi, eta)
    [~, dN, detJ] = shapeFunctionsAndJacobian(nodeCoords, xi, eta);
    B = zeros(3, 8);
    for i = 1:4
        B(1, 2*i-1) = dN(i, 1);
        B(2, 2*i) = dN(i, 2);
        B(3, 2*i-1) = dN(i, 2);
        B(3, 2*i) = dN(i, 1);
    end
end
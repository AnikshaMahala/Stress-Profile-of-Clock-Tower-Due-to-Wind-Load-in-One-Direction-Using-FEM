function [N, dNdxi, detJ] = shapeFunctionsAndJacobian(nodeCoords, xi, eta)
    N = 1/4 * [(1 - xi) * (1 - eta);
               (1 + xi) * (1 - eta);
               (1 + xi) * (1 + eta);
               (1 - xi) * (1 + eta)];

    dNdxi = 1/4 * [- (1 - eta), - (1 - xi);
                    (1 - eta), - (1 + xi);
                    (1 + eta),   (1 + xi);
                   - (1 + eta),   (1 - xi)];

    J = dNdxi' * nodeCoords;
    detJ = det(J);
end
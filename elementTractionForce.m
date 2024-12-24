function fe_traction = elementTractionForce(nodeCoords, traction, edge)
    gaussPoints = [-1 / sqrt(3), -1 / sqrt(3); 
                    1 / sqrt(3), -1 / sqrt(3); 
                   -1 / sqrt(3),  1 / sqrt(3); 
                    1 / sqrt(3),  1 / sqrt(3)];
    weights = [1, 1, 1, 1];
    fe_traction = zeros(8, 1);
    for i = 1:size(gaussPoints, 1)
        xi = gaussPoints(i, 1);
        eta = gaussPoints(i, 2);
        % 1. Compute shape functions and Jacobian at integration point
        [N, ~, detJ] = shapeFunctionsAndJacobian(nodeCoords, xi, eta);
        
        % 2. Initialize the force components for the 8x1 vector
        Fx1 = 0;  Fx2 = 0; Fx3 = 0;  Fx4 = 0;
        N(1)=0;
        N(2)=0;
        % 3. Calculate the body forces in x-direction for each node
        % Assuming bodyForce is a 2x1 matrix, i.e., bodyForce = [Fx; Fy] for each node
        % First, compute the x-components of the body force
        Fy1 = N(1) * traction * detJ*weights(i);  % Force at node 1 in x-direction
        Fy2 = N(2) * traction * detJ*weights(i);  % Force at node 2 in x-direction
        Fy3 = N(3) * traction * detJ*weights(i);  % Force at node 3 in x-direction
        Fy4 = N(4) * traction * detJ*weights(i);  % Force at node 4 in x-direction
   
        
        % 4. Update the body force vector
        fe_traction = fe_traction + [Fx1; Fy1; Fx2; Fy2; Fx3; Fy3; Fx4; Fy4];
        
        % Display the updated body force vector
        disp('Updated traction force vector:');
        disp(fe_traction);

    end
end




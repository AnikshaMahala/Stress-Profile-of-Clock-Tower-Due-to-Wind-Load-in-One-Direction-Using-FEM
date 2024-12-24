% Generalized Finite Element Analysis for Multiple Elements (Plane Stress)

% Input: Define Nodes and Element Connectivity
clc 
clear all
nodes =[
    -10.,           5.;
     5.,            5.;
     5.,    8.33333302;
   -14.2679167,    5.68240786;
     5.,           15.;
    -16.,           15.;
    -30.,           15.;
    -30.,            5.;
    -30.,           -5.;
   -16.666666,     -5.;
   -10.,           -5.;
    -5.,            5.;
     0.,            5.;
   -1.02187502,    8.01851845;
   -7.16932297,    7.43883133;
     5.,    11.666667;
    -2.,           15.;
    -9.,           15.;
   -15.3505735,    10.4836807;
   -23.,           15.;
   -30.,           10.;
   -22.4951057,     5.2550931;
   -30.,            0.;
   -23.333334,     -5.;
   -15.9750519,    0.205613613;
   -10.,            0.;
   -1.60234368,   11.4348955;
   -8.28056049,   11.0893517;
   -22.7114182,   10.1846933;
   -22.9508724,    0.115176678
];



% Element Connectivity (multiple elements)
elements = [ 1, 12, 15, 4;
12, 13, 14, 15;
13,  2,  3, 14;
3, 16, 27, 14;
16,  5, 17, 27;
14, 27, 28, 15;
27, 17, 18, 28;
15, 28, 19,  4;
28, 18,  6, 19;
6, 20, 29, 19;
20,  7, 21, 29;
19, 29, 22,  4;
29, 21,  8, 22;
8, 23, 30, 22;
23,  9, 24, 30;
22, 30, 25,  4;
30, 24, 10, 25;
10, 11, 26, 25;
25, 26,  1,  4];  % Element 3

traction_edges = containers.Map(...
    {11, 10, 9, 7, 5}, ... % Keys
    {[6, 20], [20, 7], [18, 6], [17, 18], [5, 17]}); % Values

% Input: Material Properties
E = 2e11; % Young's modulus (Pa)
nu = 0.3; % Poisson's ratio
D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2]; % Plane stress constitutive matrix

% Input: Boundary Conditions
fixedDOF = [13,14,15,16,17,18,41,42,45,46]; % 

% Input: External Forces
traction = -2e6; % Example: Traction force in y-direction along an edge
bodyForce = -9.8*7850*5; % Example: Body force (e.g., gravity) in y-direction (5 is for area)

% Step 1: Initialize Global Stiffness Matrix and Force Vector
numNodes = size(nodes, 1);
strain_at_nodes = zeros(numNodes, 3);  % For storing [epsilon_x, epsilon_y, gamma_xy] at each node
stress_at_nodes = zeros(numNodes, 3); 
numDOF = 2 * numNodes;
K = zeros(numDOF, numDOF);
F = zeros(numDOF, 1);

% Step 2: Assemble Global Stiffness Matrix and Force Vector for Multiple Elements
for elem = 1:size(elements, 1)
    elementNodes = elements(elem, :);
    nodeCoords = nodes(elementNodes, :);
    
    % Element stiffness matrix
    Ke = elementStiffness(nodeCoords, D);
    
    % Element body force (due to body forces like gravity)
    fe_body = elementBodyForce(nodeCoords, bodyForce);
    if isKey(traction_edges, elem)
    % Element traction force (traction applied on an edge)
        fe_traction = elementTractionForce(nodeCoords, traction,[1,2]); % Apply traction on edge [1, 2]
    else
        fe_traction=zeros(8, 1);
    end 
    
    % Global degrees of freedom (DOF) corresponding to element nodes
    dof = [2 * elementNodes(1) - 1, 2 * elementNodes(1), ...
           2 * elementNodes(2) - 1, 2 * elementNodes(2), ...
           2 * elementNodes(3) - 1, 2 * elementNodes(3), ...
           2 * elementNodes(4) - 1, 2 * elementNodes(4)];
       
    % Assemble the global stiffness matrix and force vector
    K(dof, dof) = K(dof, dof) + Ke;
    F(dof) = F(dof) + fe_body + fe_traction;
end


% Step 3: Apply Boundary Conditions and Solve the System of Equations
freeDOF = setdiff(1:numDOF, fixedDOF);
U = zeros(numDOF, 1); % Initialize displacement vector
U(freeDOF) = K(freeDOF, freeDOF) \ F(freeDOF); % Solve for free DOFs

fprintf('Displacement (scaled by 1e7):\n');
disp(1e2 * U); % Scale displacements for readability

% Step 4: Compute Strains and Stresses at Gauss Points for Each Element
gaussPoints = [-1 / sqrt(3), -1 / sqrt(3); 
                1 / sqrt(3), -1 / sqrt(3); 
               -1 / sqrt(3),  1 / sqrt(3); 
                1 / sqrt(3),  1 / sqrt(3)];

disp('Strains and Stresses at Gauss Points for Each Element:');
for elem = 1:size(elements, 1)
    elementNodes = elements(elem, :);
    nodeCoords = nodes(elementNodes, :);
    
    for i = 1:size(gaussPoints, 1)
        xi = gaussPoints(i, 1);
        eta = gaussPoints(i, 2);
        
        % Compute B matrix and Jacobian determinant
        [B, detJ] = strainDisplacementMatrix(nodeCoords, xi, eta);
        
        % Extract element displacement vector (for the current element)
        dof = [2 * elementNodes(1) - 1, 2 * elementNodes(1), ...
               2 * elementNodes(2) - 1, 2 * elementNodes(2), ...
               2 * elementNodes(3) - 1, 2 * elementNodes(3), ...
               2 * elementNodes(4) - 1, 2 * elementNodes(4)];
        Ue = U(dof);
    
        % Compute strain and stress for the current Gauss point
        strain = B * Ue;
        stress = D * strain;
        for j = 1:length(elementNodes)
            nodeIndex = elementNodes(j);  % Get the node index
            
            % Accumulate strain and stress components for the node
            strain_at_nodes(nodeIndex, :) = strain_at_nodes(nodeIndex, :) + strain';
            stress_at_nodes(nodeIndex, :) = stress_at_nodes(nodeIndex, :) + stress';
        end
        
        % Display results for each element and Gauss point
        fprintf('Element %d, Gauss Point (xi = %.2f, eta = %.2f):\n', elem, xi, eta);
        fprintf('Strain (scaled by 1e7):\n');
        disp(1e2* strain);
        fprintf('Stress (scaled by 1e7):\n');
        disp(1e2 * stress);
    end
end

% Define scaling factors for displacements
scaleFactor = 4000; % Adjust for better visualization

% Compute new displacement UX and UY
deltaUx = scaleFactor * U(1:2:end); % Extract x-components of displacement
deltaUy = scaleFactor * U(2:2:end); % Extract y-components of displacement
Ux = nodes(:, 1); % Original X-coordinates of nodes
Uy = nodes(:, 2); % Original Y-coordinates of nodes

% Updated displacements
UX = Ux + deltaUx;
UY = Uy + deltaUy;

% Updated node coordinates after displacement
deformedNodes = [UX, UY];

% Plot original mesh
figure;
hold on;
axis equal;
grid on;
title('Original and Deflected Mesh');
xlabel('X-coordinate');
ylabel('Y-coordinate');

% Plot the undeformed (original) mesh
for elem = 1:size(elements, 1)
    elementNodes = elements(elem, :);
    originalCoords = nodes(elementNodes, :);
    h1 = plot([originalCoords(:, 1); originalCoords(1, 1)], ...
              [originalCoords(:, 2); originalCoords(1, 2)], ...
              'b-', 'LineWidth', 1.5); % Blue lines for original mesh
end

% Plot the deformed (deflected) mesh
for elem = 1:size(elements, 1)
    elementNodes = elements(elem, :);
    deformedCoords = deformedNodes(elementNodes, :);
    h2 = plot([deformedCoords(:, 1); deformedCoords(1, 1)], ...
              [deformedCoords(:, 2); deformedCoords(1, 2)], ...
              'r-', 'LineWidth', 1.5); % Red lines for deformed mesh
end
% Annotate element numbers
for elem = 1:size(elements, 1)
    elementNodes = elements(elem, :);
    elementCoords = nodes(elementNodes, :);
    centroid = mean(elementCoords, 1); % Compute centroid of the element
    text(centroid(1), centroid(2), sprintf('%d', elem), ...
         'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Add legend with specific handles
legend([h1, h2], {'Original Mesh', 'Deflected Mesh'}, 'Location', 'Best');



hold off;

% Average strain and stress values over each node (optional)
strain_at_nodes = strain_at_nodes / 4;  % Divide by the number of Gauss points (4 in this case)
stress_at_nodes = stress_at_nodes / 4;  % Divide by the number of Gauss points (4 in this case)

% Plot Stress components
figure;
subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), stress_at_nodes(:, 1), 'EdgeColor', 'none');
colorbar;
title('Stress - \sigma_x');
xlabel('X');
ylabel('Y');
zlabel('Stress');
view(2);  % View from above
axis equal;
shading interp;

subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), stress_at_nodes(:, 2), 'EdgeColor', 'none');
colorbar;
title('Stress - \sigma_y');
xlabel('X');
ylabel('Y');
zlabel('Stress');
view(2);
axis equal;
shading interp;

subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), stress_at_nodes(:, 3), 'EdgeColor', 'none');
colorbar;
title('Stress - \tau_{xy}');
xlabel('X');
ylabel('Y');
zlabel('Stress');
view(2);
axis equal;
shading interp;

% Plot Strain components
figure;
subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), strain_at_nodes(:, 1), 'EdgeColor', 'none');
colorbar;
title('Strain - \epsilon_x');
xlabel('X');
ylabel('Y');
zlabel('Strain');
view(2);
axis equal;
shading interp;

subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), strain_at_nodes(:, 2), 'EdgeColor', 'none');
colorbar;
title('Strain - \epsilon_y');
xlabel('X');
ylabel('Y');
zlabel('Strain');
view(2);
axis equal;
shading interp;

subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), strain_at_nodes(:, 3), 'EdgeColor', 'none');
colorbar;
title('Strain - \gamma_{xy}');
xlabel('X');
ylabel('Y');
zlabel('Strain');
view(2);
axis equal;
shading interp;



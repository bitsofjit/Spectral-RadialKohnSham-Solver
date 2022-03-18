function [Vh] = solveradialpoisson(xe, connect, numel, ...
                                   nodaldofs, dofinfo, numdofs, ...
                                   xq, wtq, phiq, dphiq, rhoq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Inputs
% ====== 
% rhoq(:,)               : quadrature-grid representatin of electron density; 
%                          rho(i,j) = rho value at i-th quadrature point 
%                          of j-th element
% phiq(:,:), dphiq(:,:)  : parent basis fn values at 
%                          quadrature points xiq; 
%                          phihq(i,j) = value of j-th func
%                          at i-th quadrature point
% 
% Outputs
% =======  
%

m = 0;                      % size of index arrays

for e = 1:numel
    conn = connect(e,:); 
    sctr = [ nodaldofs(conn).list ];
    m = m + (length(sctr))^2;
end

I = zeros(m,1); 
J = I; 
XK = I; 
f = zeros(numdofs, 1);
%
% compute and assemble Stiffness matrix K and load vector f
% 
ind = 1; % indices for arrays I, J, and XK

%
% loop over all the elements
%
for e = 1:numel
    
    conn = connect(e,:);
    sctr = [nodaldofs(conn).list];
    
    dofs_to_nodes = dofinfo(sctr);      % dofs = DOFs of the element
                                        % global node Ids 
    
    detj = (xe(e+1) - xe(e))/2.0;
    %
    % compute element contributions to Stiffness matrix and load vector:
    %
    [ke, fe] = radialpoissonelemeqns(conn, dofs_to_nodes, phiq, dphiq, ...
                                     rhoq(:,e), xq(:,e), wtq, detj);
    
    %
    % insert ke, fe into vector arrays
    %
    sctrlen = length(sctr);
    temp = meshgrid(sctr);
    I(ind:ind+sctrlen^2-1) = temp(:);
    temp = temp';
    J(ind:ind+sctrlen^2-1) = temp(:);
    
    XK(ind:ind+sctrlen^2-1) = ke(:);
    f(sctr) = f(sctr) + fe;
    ind = ind + sctrlen^2;
    
end

%
% Create Stiffness matrix from the array
%

fprintf('Performing sparse assembly of Stiffness Matrix (%d x %d) . . .\n\n',numdofs,numdofs);

bigk = sparse(I, J, XK, numdofs, numdofs);

clear I; clear J; clear XK;


%
% solve the linear system Ku = f
%
fprintf('Solving the linear system [K]{u} = {b} . . .\n\n');
Vh = bigk\f;

return
end

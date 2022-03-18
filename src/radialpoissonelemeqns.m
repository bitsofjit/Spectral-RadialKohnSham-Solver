function [ke, fe] = radialpoissonelemeqns(nodes, dofs_to_nodes, ...
                                          phiq, dphiq, rhoq, xqe, wtq, detj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
%
% Inputs
% ======                               
% rhoq(:)          : quadrature point values of the electron density.
%                    rhoq(i) is the value of density at the i-th quadrature 
%                    point of the element under consideration
%
% Outputs
% =======  
%
ndofs = length(dofs_to_nodes);
dofs_to_locindices = zeros(1,ndofs);

for i = 1:ndofs
  dofs_to_locindices(i) = find(nodes == dofs_to_nodes(i));
end


%
% initialize the stiffness matrix and load vector to zero
%
ke = zeros(ndofs,ndofs);
fe = zeros(ndofs,1);

%
% loop over Gauss points in xi-direction
%

for i = 1:length(xqe)
    r = xqe(i);
    %
    % FE shape functions and local derivatives
    %
    NFE     = phiq(i,:);
    dNFEdxi = dphiq(i,:);
    
    N = NFE(dofs_to_locindices);
    dNdxi = dNFEdxi(dofs_to_locindices);
    
    % derivative of shape function with respect to global coordinate r
    dNdr = dNdxi/detj;
    
    %
    % The element Stiffness matrix (ke) and load vector (fe):
    %    
    % Considering Vp = r * Vh as the varibale to solve:
    ke = ke + detj*wtq(i)*(dNdr'*dNdr);
    fe = fe + detj*wtq(i)*N'*r*rhoq(i);
end

return
end

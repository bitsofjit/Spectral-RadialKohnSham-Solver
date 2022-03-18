function [he, se] = radialschroedingerelemeqns(nodes, dofs_to_nodes, ...
                                               Veffelem, l, phiq, dphiq, ...
                                               xqe, wtq, detj)
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
% Veff(:)          : VH[rho(:)] + Vxc[rho(:)] + Vext(:);
%                    A vector of length nsp. The i-th elements of Veff
%                    vector contains the sum of the values of VH, 
%                    Vxc, and Vext evaluated at  the i-th quadrature 
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
% initialize element Hamiltonian (he) and overlap (se) matrices to zero
%
he = zeros(ndofs,ndofs);
se = zeros(ndofs,ndofs);

%
% loop over Gauss points in xi-direction
%
for i = 1:length(xqe)
    
    r = xqe(i);    
    %
    % FE shape functions and local derivatives
    %
    NFE = phiq(i,:);
    dNFEdxi = dphiq(i,:);
    
    N = NFE(dofs_to_locindices);
    dNdxi = dNFEdxi(dofs_to_locindices);
    
    % derivative of shape function with respect to global coordinate r
    dNdr = dNdxi/detj;
    
    %
    % The element Hamiltonian (he) and overlap matrices (se)
    %
    mass = N'*N;
    
    % Using u = P_{nl} = r * R_{nl} as the varibale to solve:
    he = he + detj*wtq(i)*(0.5*(dNdr'*dNdr) + ...
                          (Veffelem(i) + l*(l + 1)/2.0/r/r)*mass);
    
    se = se + detj*wtq(i)*mass;
end

return
end

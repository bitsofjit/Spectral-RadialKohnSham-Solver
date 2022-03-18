function uq = fetoquad(xi, xiq, connect, ufe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% transforms ufe from FE-coefficient to quadrature-grid representation
% ufe is a full FE coefficient vector, having values for all nodes in 
% the mesh, including domain-boundary nodes
% 
% Inputs
% ======
% ufe(:)           : full FE coefficient vector, having values for all nodes 
%                    in a mesh, including domain-boundary nodes.
% xi(:)            : parent basis nodes; xi(i) = coordinate of i-th 
%                    parent basis node
% xiq(:)           : quadrature points
% connect(:,:)     : nodal connectivity; connect(i,j) = index of basis node
%                    corresponding to local node j of element i
% 
% Output
% ======        
% uq(:,:)          : quadrature-grid representatin of ufe;
%                    uq(i,j) = value at i-th quadrature point of j-th element

% Tabulate parent basis at quadrature points:
nparntnodes = length(xi);
p = nparntnodes - 1;
nsp = length(xiq);
numel = size(connect,1);

[phiq, ~] = getshapefunc(p, xiq);

%  Evaluate at quad points in each element
uq = zeros(nsp,numel);
for e = 1:numel
    ue = ufe(connect(e,:));
    for iq = 1:nsp
        uq(iq, e) = dot(phiq(iq,:), ue);
    end
end

end
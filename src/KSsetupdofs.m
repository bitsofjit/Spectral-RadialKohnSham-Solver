function [nodaldofs, dofinfo, numdofs] = KSsetupdofs(numnod)

isboundarynode = zeros(1,numnod);

% Dirichlet BCs on boundary nodes
isboundarynode(1) = 1; 
isboundarynode(numnod) = 1;

nodaldofs(numnod).list = [];
numdofs = 0;

for n = 1:numnod
  if (~isboundarynode(n))
    numdofs = numdofs + 1;            % classical DOF
    nodaldofs(n).list = [numdofs];
    dofinfo(numdofs) = n;             % node number for dof equal to numdofs
  end

end

return
end

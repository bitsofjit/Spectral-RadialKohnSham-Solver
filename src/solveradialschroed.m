function [lambda, eigv, eigsflag] = solveradialschroed(xe, connect, numel, ...
                            nodaldofs, dofinfo, numdofs, ...
                            Veffq, l, xq, xiq, wtq, phiq, dphiq, nbeigs, eigopts)
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Inputs
% ====== 
% Veffq(:,:)     : quadrature-grid representatin of Veff - 
%                  Veffq(:,:) := VH[rhoq(:,:)] + Vxc[rhoq(:,:)] + Vext(:,:);
%                  Veffq(:,:) = value at i-th quadrature point of j-th element
% 
% phiq(:,:), dphiq(:,:)  : parent basis fn values at 
%                          quadrature points xiq; 
%                          phihq(i,j) = value of j-th func
%                          at i-th quadrature point
% xq(:,:)                : quadrature pts in physical space; 
%                          xq(i,j) = coordinate of i-th point in j-th element
% Outputs
% =======  
%
% 
global quadtype;
m = 0; % size of index arrays

for e = 1:numel
  conn = connect(e,:); sctr = [ nodaldofs(conn).list ];
  m = m + (length(sctr))^2;
end

I = zeros(m,1); J = I; XH = I; XS = I;

%
% compute and assemble Hamiltonian matrix H and overlap matrix S
%
% loop over all the elements
%

ind = 1; % indices for arrays I, J, XH, XS

for e = 1:numel
    
    conn = connect(e,:);
    sctr = [nodaldofs(conn).list];
    
    dofs_to_nodes = dofinfo(sctr);      % dofs = DOFs of the element
                                        % global node Ids
    detj = (xe(e+1) - xe(e))/2.0;
    %
    % compute element contributions to Stiffness matrix and load vector:
    %                                    
    [he, se] = radialschroedingerelemeqns(conn, dofs_to_nodes, ...
                                          Veffq(:,e), l, phiq, dphiq, ...
                                          xq(:,e), wtq, detj);
    %
    % insert he, se into vector arrays
    %
    sctrlen = length(sctr);
    
    temp = meshgrid(sctr);
    I(ind:ind+sctrlen^2-1) = temp(:);
    temp = temp';
    J(ind:ind+sctrlen^2-1) = temp(:);
    XH(ind:ind+sctrlen^2-1) = he(:);
    XS(ind:ind+sctrlen^2-1) = se(:);
    ind = ind + sctrlen^2;
    
end

%
% Create H and S-matrices from the arrays
%

fprintf('Performing sparse assembly of system matrices (%d x %d) . . .\n\n',numdofs,numdofs);

bigh = sparse(I, J, XH, numdofs, numdofs);
bigs = sparse(I, J, XS, numdofs, numdofs);

clear I; clear J; clear XH; clear XS;

%
% solve the eigen-problem (generalized/standard): H c = \lambda S c
%
% Symmetrize (in case bigs is nearly-symmetric due to round-off, etc.)
bigs = 0.5*(bigs + bigs');

% opt.TOLERANCE    = 1.0E-10;
% opt.MAXITERATION = 1000;
% opt.iluThresh    = 0.1; 
% 
% [lambda, eigv] = eigifp(bigh, bigs, nbeigs, opt);

if (size(phiq,2) == length(xiq)) && (quadtype == 2)    % nsp = p + 1
    % Call eigs as a standard eigenproblem solver:
    fprintf('Solving standard eigenproblem: [H]{c} = e{c} . . .\n\n');
    [eigv, D, eigsflag] = eigs(bigs\bigh, nbeigs , 'SM', eigopts);
else
    % Call eigs as a generalized eigenproblem solver:
    fprintf('Solving generalized eigenproblem: [H]{c} = e[S]{c} . . .\n\n');
    [eigv, D, eigsflag] = eigs(bigh, bigs, nbeigs , 'SA', eigopts);
%     DON't do this as well
%     [eigv, D, eigsflag] = eigs(bigs\bigh, nbeigs , 'SA', eigopts);   %


%   WARNING: Setting sigma = 0 screws up the whole problem ****
%   DO NOT USE THAT!
%   [eigv, D, eigsflag] = eigs(bigh, bigs, nbeigs , 0.0, eigopts);

end

if eigsflag ~= 0
    fprintf('\n')
    fprintf('       **************************      \n')
    fprintf('             FATAL ERROR! \n')
    fprintf('       **************************      \n')
    fprintf('\n')
    fprintf('All the required eigenvalues did not converge! \n');
    fprintf('\n')
    fprintf('Consider changing eigopts.tol/eigopts.maxit and re-run. \n')
    condbigs = condest(bigs);
    condbigh = condest(bigh);
    fprintf('\n')
    fprintf('The Condition number of Overlap matrix at exit %-16.10e \n', condbigs) 
    fprintf('\n')
    fprintf('and the Condition number of Hamiltonian matrix at exit %-16.10e \n', condbigh)
    fprintf('\n')
    fprintf('Aborting . . . \n');
    fprintf('\n')
    lambda = NaN;
    return
else
    fprintf('\n')
    fprintf('All the required eigenvalues converged within set tolerance. \n');
    fprintf('\n')
    lambda = diag(D);
    return
end



end
                                    
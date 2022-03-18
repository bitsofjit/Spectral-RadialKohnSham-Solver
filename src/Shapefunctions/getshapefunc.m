function [phiq, dphiq] = getshapefunc(p, xiq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Tabulate parent basis at quadrature points:
% 
% Inputs
% ====== 
% 
% 
% Outputs
% =======  
nsp = length(xiq);

phiq = zeros(nsp,p+1);    % parent basis fn values at 
%                           quadrature points xiq; 
%                           phihq(i,j) = value of j-th func
%                           at i-th quadrature point
dphiq = zeros(nsp,p+1);   % derivative of the parent basis fn values at 
%                           quadrature points xiq; 
%                           dphihq(i,j) = value of j-th func
%                           at i-th quadrature point

for iq = 1:nsp
    [phiq(iq,:), dphiq(iq,:)] = shapefunctions(p, xiq(iq)); % all the shapefunctions
end
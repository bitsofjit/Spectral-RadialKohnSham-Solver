%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [NFE dNdxi] = shapefunctions(eorder,xi) 
% Purpose
% =======
% Compute the shape function and its local derivative at the Gauss point
%

function [NFE, dNdxi] = shapefunctions(eorder,xi) 
%
% shape function
%
% nodes 1 to eorder+1
%

n = eorder+1;

[NFE, dNdxi] = lagrange(xi,n);

return
end

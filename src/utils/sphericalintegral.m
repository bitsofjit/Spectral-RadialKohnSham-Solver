function r = sphericalintegral(xe, xq, wtq, fq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Spherical integral of function f(r) (no \theta, \phi dependence) 
% defined by values [fq] at the quadrature points on mesh xe. 
% i.e., \int_0^rmax 4 pi r^2 f(r) dr, where rmax = max r value in mesh


r = integrateonmesh(xe, wtq, 4.0.*pi.*xq.^2.*fq);

return
end

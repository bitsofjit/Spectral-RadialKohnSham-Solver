function [T_s, E_ee, E_en, E_xc, E_tot] = computetotalenergy(xe, xq, wtq, ...
                                          rhooutq, Veffq, Vhq, excq, ...
                                          Z, E_band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% =======      
% Variational, quadratically convergent form of total energy
% Returns the total energy as well as its parts
%
% Inputs
% ======
% xe(:)          : xe(i/i+1) = coord of left/right boundary of ith element
% xq(:,:)        : quadrature pts in physical space; 
%                  xq(i,j) = coordinate of i-th point in j-th element
% wtq(:)         : Quadrature weights
% rhooutq(:,:)   : radial density at quadrature points; rhoq(i,j) = value
%                  at i-th quadrature point of j-th element
% veffq(:, :)    :
% V_h_q(:, :)    :
% e_xc_q(:, :)   :
% E_band         :
% Z              :
% 
% Output
% ======
% 

xqsqrd = xq.^2;

T_s = E_band + 4.0*pi*integrateonmesh(xe, wtq, Veffq.*rhooutq.*xqsqrd);

E_ee = -2.0*pi*integrateonmesh(xe, wtq, Vhq.*rhooutq.*xqsqrd);

E_en = 4.0*pi*Z*integrateonmesh(xe, wtq, rhooutq.*xq);

E_xc = -4.0*pi*integrateonmesh(xe, wtq, excq.*rhooutq.*xqsqrd);

E_tot = T_s + E_ee + E_en + E_xc;

return

end
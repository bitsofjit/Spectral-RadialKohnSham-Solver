function rho = thomasfermidensity(r, Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Generalized Thomas-Fermi density (charge)
% Inputs
% ====== 
% r(:)      : Radial grid
% Z         : Atomic number
% 
% Output
% ======
% rho(:)    : Thomas-Fermi charge density evaluated at the input grid points

V = thomasfermipotential(r, Z, false);

rho = -1.0./3.0/pi^2.*(-2.0*V).^(3.0/2.0);      % Charge density

return
                                                                                                                                                                     
end

function V = thomasfermipotential(r, Z, cut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Generalized Thomas-Fermi atomic potential
% 
% Inputs
% ====== 
% r(:)      : Radial grid
% Z         : Atomic number
% cut       : Cut the potential
% 
% Output
% ======
% V(:)      : Thomas-Fermi potential evaluated at the input grid points
% 
if nargin < 3
    cut = true;
end

x = r.*(128.0*Z/9.0/pi^2).^(1.0/3.0);

%  Z_eff(x) = Z * phi(x), where phi(x) satisfies the Thomas-Fermi equation:
%                   phi'' = phi^(3/2) / sqrt(x)
%  with boundary conditions:
%    phi(0)  = 1
%    phi(oo) = 0
%  There is no analytic solution, but one can solve this approximately. We use:
%  http://arxiv.org/abs/physics/0511017

% Parameters:
alpha =  0.7280642371;
beta  = -0.5430794693;
gamma =  0.3612163121;

sqrtx = sqrt(x);
Zeff = Z.*(1 + alpha*sqrtx + beta*x.*exp(-gamma*sqrtx)).^2.* ...
           exp(-2.0*alpha*sqrtx);
    
%  This keeps all the eigenvalues of the radial problem negative:
if cut == 1
    Zeff(Zeff < 1.0) = 1.0;
end
    
V = -Zeff./r;

return
end

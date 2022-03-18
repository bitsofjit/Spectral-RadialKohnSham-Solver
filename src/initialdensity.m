function rhoq = initialdensity(Q, xe, xq, wtq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Constructs localized charge density of total charge -Q at quadrature points
% 
% Inputs
% ====== 
% Q         : total charge; density of total charge -Q is constructed
% xe(:)     : xe(i/i+1) = coord of left/right boundary of ith element
% xiq(:)    : quadrature points
% wtq(:)    : quadrature weights
% xq(:,:)   : quadrature pts in physical space; 
%             xq(i,j) = coordinate of i-th point in j-th element
% 
% Output
% ======
% rhoq(:,:) : radial charge density at quadrature points; 
%             rhoq(i,j) = value at i-th quadrature point of j-th element
% 
numel = length(xe) - 1;

rhoq = zeros(size(xq));

% rhoq = zeros(length(xi),numel);

for e = 1:numel          % tabulate over elements
    xqe = xq(:,e);
    % rhoq(:,e) = exp(-xqe.^2./4.0)             % Gaussian
    % rhoq(:,e) = 1.0./cosh(xqe).^2             % 1/cosh^2
    %
    % - ( Thomas-Fermi charge-density ) ==> electronic-density
    %
    rhoq(:,e) = -thomasfermidensity(xqe, Q);
end

integ = sphericalintegral(xe, xq, wtq, rhoq);
rhoq = -rhoq./integ*Q;

return

end

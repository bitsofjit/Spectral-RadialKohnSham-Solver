function sm = integrateonmesh(xe, wtq, fq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Integral of function f(x) defined by values fq at quadrature points 
% on mesh xe, without additional factors; i.e., \int_xmin^xmax f(x) dx
% 
% Inputs
% ====== 
% xe(:)     : xe(i/i+1) = coord of left/right boundary of i-th element
% wtq(:)    : quadrature weights
% fq(:,:)   : integrand values at quadrature points; 
%             fq(i,j) = func. value at i-th quad pt in j-th element
% 
% Output
% ======
% 
numel = length(xe) - 1;
sm = 0;

for e = 1:numel            % sum over elements
    % Jacobian of affine transformation from parent coords
    % xi in [-1,1] to coords x in element [xa,xb]: 
    % x = al xi + be, al=(xb-xa)/2, be=(xb+xa)/2     
    jac = 0.5*abs(xe(e+1) - xe(e));
    sm = sm + sum(wtq'.*fq(:,e))*jac;   % dot(wtq, fq(:,e)*jac);
end

return
end
function xq = getphysicalquadpts(xe, xiq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Generates quadrature points xq in physical coordinate system via 
% affine mapping of parent quad points xiq to elements xe
% 
% Inputs
% ====== 


% xe(:)     : xe(i/i+1) = coord of left/right boundary of ith element
% xiq(:)    : quadrature points in [-1, 1] (parent domain)
% 
% Output
% ======
% xq(:,:)   : quadrature pts in physical space; 
%             xq(i,j) = coordinate of i-th point in j-th element


numel = length(xe) - 1;     nsp = length(xiq);

xq = zeros(nsp,numel);    % assumes same nsp for all the elements

for e = 1:numel
   xl = xe(e); xr = xe(e+1);
   % affine transformation: [-1 1] --> [xl xr]
   xq(:,e) = (xr - xl)/2.0*xiq + (xr + xl)/2.0;
end


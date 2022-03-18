function [excq, Vxcq] = exhangecorrelation(type, rhoq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% Calculates exchange-correlation (LDA density and potential so far) 
% from the electron density rhoq. rhoq must be positive
% 
%
% Inputs
% ======  
% rhoq(:,:)    : quadrature-grid representatin of rho; 
%                rho(i,j) = rho value at i-th quadrature point 
%                of j-th element
% Output
% ======        
% Vxcq(:,:)    : quadrature-grid representatin of Vxc;
%                Vxc(i,j) = value at i-th quadrature point of j-th element

type = lower(type);     % LDA density and potential so far

switch type
    case 'lda'          % The Local Density Approximation (LDA)
        % Parameters:
        y0 = -0.10498;
        b  = 3.72744;
        c  = 12.9352;
        A  = 0.0621814;
        
        Q  = sqrt(4.0*c - b^2);
        rs = (3.0./(4.0*pi*rhoq)).^(1.0/3.0);
        y  = sqrt(rs);
        Y  = computeY(y, b, c);
        Y0 = computeY(y0, b, c);
        
        % Vosko-Wilk-Nusair correlation term:
        ec = A/2.0*(log(y.^2./Y) + 2.0*b/Q.*atan(Q./(2.0*y + b)) - ...
             b*y0./Y0.*(log((y - y0).^2./Y) + 2.0*(b + 2.0*y0)./ ...
             Q.*atan(Q./(2.0*y + b))));
         
        Vc = ec - A/6.0*(c*(y - y0) - b*y0*y)./((y - y0).*Y);
        
        % Electron gas exchange term:
        ex = -3.0./(4.0*pi)*(3.0*pi^2*rhoq).^(1.0/3.0);
        Vx = 4.0*ex./3.0;
        
        excq = ex + ec;
        Vxcq = Vx + Vc;
        
    otherwise
        error('Wrong exchange-correlation type entered. Use lda instead');
end

rhonegval = any(rhoq(:) < 0.0);



if rhonegval
    disp('Density is negative!');
    indcs = find(rhoq(:) < 0.0);
    excq(indcs) = 0.0;
    Vxcq(indcs) = 0.0;
end

return

end


function Y = computeY(y, b, c)
Y = y.^2 + b.*y + c;
end


        
        
        
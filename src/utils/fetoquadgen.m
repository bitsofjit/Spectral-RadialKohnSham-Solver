function uqII = fetoquadgen(xiI, xiqII, xeI, xeII, connectI, ufeI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% ======= 
% transforms ufe from FE-coefficient of mesh-I to quadrature-grid 
% representation of mesh-II. 
% ufe is a full FE coefficient vector, having values for all nodes in 
% mesh-I, including the domain-boundary nodes
% 
% Inputs
% ======
% ufeI(:)          : full FE coefficient vector, having values for all nodes 
%                    in a mesh-I, including domain-boundary nodes.
% xiI(:)           : parent basis nodes; xi(i) = coordinate of i-th 
%                    parent basis node of mesh-I
% xiqII(:)         : quadrature points of mesh-II
% connect(:,:)     : nodal connectivity; connect(i,j) = index of basis node
%                    corresponding to local node j of element i
% xq(:,:)          : quadrature pts in physical space; 
%                    xq(i,j) = coordinate of i-th point in j-th element
% 
% Output
% ======        
% uqII(:,:)          : quadrature-grid representatin of ufe on mesh-II;
%                      uqII(i,j) = value at i-th quadrature point of 
%                      j-th element

% Generate global coordinates of quadrature points of mesh-II:
nspII = length(xiqII);          % # of quadrature points of mesh-II 
xqII = getphysicalquadpts(xeII, xiqII);

numelI  = length(xeI) - 1;      % # of elements of mesh-I
numelII = length(xeII) - 1;     % # of elements of mesh-II
pI = length(xiI) - 1;            % element order of mesh-I

% Evaluate at quad points of mesh-II in each element
% Initialize output
uqII = zeros(nspII,numelII);

for eII = 1:numelII              % loop over elements of mesh-II
    % Physical coordintaes of quadrature-points inside element eII in mesh-II
    xqIIeII = xqII(:,eII);       % length nspII
    
    for iqII = 1:length(xqIIeII) % loop over q-points of eII
        xsampleII = xqIIeII(iqII);
        
        % Now loop over elements of mesh-I to scan which element (can be
        % only one) contains quadrature point xsampleII
        for eI = 1:numelI        
            if iswithin(xsampleII, xeI(eI), xeI(eI+1)) == 1
                % get the FE-nodal solution for that element
                ueI = ufeI(connectI(eI,:));
                
                xl = xeI(eI); 
                xr = xeI(eI+1);
                l = xr - xl;
                Nl = (xr - xsampleII)/l;
                Nr = (xsampleII - xl)/l;
                % affine transformation: [xl xr] --> [-1 1]
                xisample = Nr - Nl;     % Nl*(-1) + Nr*(+1)
                
                % Tabulate parent basis at quadrature point xsampleII:
                [phi, ~] = getshapefunc(pI, xisample);
                
                % Compute the FE solution at quadrature point xsampleII
                uqII(iqII,eII) = dot(phi, ueI);
                
                % since xsampleII can be contained by just one element. 
                % Once found, hence just "break" the loop to move on 
                % to the next quadrature point
                break
            else
                continue
            end
        end
    end
end

end



function r = iswithin(x, low, hi)
% returns logical true if low < x < hi
r = (x >= low) & (x <= hi);
end

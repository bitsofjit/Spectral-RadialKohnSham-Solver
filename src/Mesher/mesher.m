%----------------------------------------------------------------------------
% Function Syntax
% ===============
% function [coord,connect] = mesh(meshtype, param, options)
%
% Inputs
% ======   
% mesh generated on [rin rout]
% mesh type: 1 = uniform, 2 = logarithmic, 3 = exponential
% 
% options             : data structure with element information 
%
% Outputs
% =======        
% coord(:)            : 1D array which stores the global coordinates of
%                       the nodes; coord(I,1) : node I coordinates
% connect(:,:)        : nodal connectivity; connect(i,j) = index of basis
%                       node corresponding to local node j of element i
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function meshsetup = mesher(meshtype, param, meshopts)

element_type = meshopts.element_type;

switch element_type
  case 'spectral'
  otherwise
    fprintf('Mesh generator not coded for element type = %s \n',element_type);
    error('Aborting . .');
end

p      = meshopts.element_order;
numel  = meshopts.numel;
domain = meshopts.domain;
plot   = meshopts.verbose;
rmin   = domain(1);
rmax = domain(2);

xin = linspace(-1.0, 1.0, p+1);
numnod = p*numel + 1;

% Checks:
if (rmax <= rmin)
    error("Error! rmax <= rmin. Aborting . ."); 
end
if (numel < 1)
    error("Error! Number of elements < 1. Aborting . .");
end

switch meshtype
    case 1      % uniform
        xe = linspace(rmin, rmax, numel+1);
        
    case 2      % logarithmic
        % using DOI: 10.1016/j.cpc.2013.02.014
        % dx_Ne/dx_1 = mpars(1) = a;
        % dx_{i+1}/dx_i = q, where dx_i = length of i-th element
        % => r = mpars(1)^(1/(Ne - 1))
        % x_i = (b - a)(q^i - 1)/(q^Ne - 1) + a, 
        % where x_i = right boundary of i-th element
        
        % Some Checks:
        if (numel < 2), error("Error! Logarithmic mesh requires numel >= 2. Aborting . . ."); end
        
        a = param;
        if (a == 1), error("Error! Logarithmic mesh requires a != 1. Aborting . . ."); end
            
        % xe(i/i+1) = coord of left/right boundary of ith element
        xe = zeros(numel+1, 1);
        xe(1,1) = rmin; xe(end,1) = rmax;
        
        q = a^(1/(numel - 1));
        
        for i = 1:numel
            xe(i+1,1) = (q^i - 1.0)/(q^numel - 1.0)*(rmax - rmin) + rmin;
        end
        disp(' ')     
    case 3      % exponential mesh (Not tested)      
        a = param;
        
        if (a < 0)
            error("Error! a > 0 required. Aborting . .");
        elseif (a == 1)
            alpha = (rmax - rmin)/numel;
            
            for i = 1:numel+1
                xe(i) = alpha*(i - 1.0) + rmin;
            end
        else    
            if (numel > 1)
                xe = zeros(numel+1,1);
                beta = log(a)/(numel - 1);
                alpha = (rmax - rmin)/(exp(beta*numel) - 1.0);
                for i = 1:numel+1 
                    xe(i) = alpha*(exp(beta*(i - 1)) - 1.0) + rmin;
                end 
            elseif (numel == 1)
                xe(1) = rmin;
                xe(2) = rmax;
            else
                error("Error! numel >= 1 required. Aborting . .");
            end
        end
              
    otherwise 
        error('More coming soon!');
end
% 
% Generate nodal coordinates in equi-spaced settings
% 
coord  = zeros(numnod, 1);
for e = 1:numel
    coord((e-1)*(p+1)+1:e*(p+1),1) = linspace(xe(e), xe(e+1), p+1);
end

%
% connectivity
%
connect = zeros(numel,p+1);
for e = 1:numel
  is = p*(e-1) + 1; 
  connect(e,:) = is:is+p;
end

%
% modify nodal coordinates (Lobatto points)
%

if (p >= 2)
  xin = lobatto_points(p+1);
  for e = 1:numel
    is = p*(e-1) + 1; 
    a = xe(e); b = xe(e+1);
%     a = coord(is,:); b = coord(is+p,:);
    for n = is:is+p
      index = n-is+1;
      coord(n,:) = 0.5*(b + a) + 0.5*(b - a).*xin(index);
    end
  end
end

meshsetup = struct('elemendcoord',xe, 'parentnodecoord',xin, 'coord',coord, 'connect',connect);


%-------------------------------------------------------------------
% Mesh plotting
%-------------------------------------------------------------------

if (strcmp(plot,'on'))
    if numel < 10 && p <= 5
        plotnodes(xe, coord);
        hold off;
    end
else
    fprintf('\n');
    fprintf('Disabling mesh-plotting because of too many elements/basis functions. \n');
    fprintf('\n');
end

return

end

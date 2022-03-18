%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [phi dphi] = lagrange(xi,n)
%
% Purpose
% =======
% Compute the Lagrange shape functions and its derivatives for xi \in [-1,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi, dphi] = lagrange(xi,n)

% n = number of shape functions (n = order+1)

xip = lobatto_points(n);

xivec = xi*ones(1,n);
vec = xivec - xip;

phi  = zeros(1,n);
dphi = zeros(1,n);
denom = zeros(1,n);

for i = 1:n
  num = 1;
  den = 1.;
  for j = 1:n
    if (j ~= i)
      num = num*vec(j);
      den = den*(xip(i) - xip(j)); 
    end
  end
  denom(i) = den;
  phi(i) = num/den;
end

for i = 1:n
  avec = vec; avec(i) = 1.;
  num = 0.;
  for j = 1:n
    if (j ~= i)
      bvec = avec; bvec(j) = 1.;
      num = num + prod(bvec);
    end
  end
  dphi(i) = num/denom(i);

end

return
end

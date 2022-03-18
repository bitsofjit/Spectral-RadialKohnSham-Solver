function plotshapefunctions(domain,m,order)
%
% Usage: plotshapefunctions(domain,m,order)
% Example: plotshapefunctions([0 1],100,3) => [0,1] is the domain
%                                             output at 100 points 
%                                             cubic shape functions
%

close all;
if (nargin < 3)
  help(mfilename);
  domain = [0 1]; m = 40; order = 3;
end

a = domain(1); b = domain(2);
n = order + 1; % no of shapefunctions

xi = linspace(-1,1,m)';
x = 0.5*(1.-xi)*a + 0.5*(1+xi)*b;
for i = 1:length(x)
  xicoord = xi(i);
  [p dp] = lagrange(xicoord,n); % p and dp are row vectors
  phi(i,:) = p'; % column
end

figure;
set(0,'DefaultAxesColorOrder',[0 0 0]);

setplotcommand(x,phi);

xlabel('x','FontSize',22);
ylabel('\phi','FontSize',22);
set(gca,'FontSize',22);

[m n ] = size(phi);
setlegend(n,'basis');

hold on;
xi_lob = lobatto_points(n);
xI = 0.5*(1.-xi_lob)*a + 0.5*(1+xi_lob)*b;
plot(xI,zeros(n,1),'mo','MarkerEdgeColor','k','MarkerFaceColor',[0.1 1 0.1],'MarkerSize',12,'linewidth',1);

return
end

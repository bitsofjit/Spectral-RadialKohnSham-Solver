%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Subhajit Banerjee & N. Sukumar
% July 2017
% UC Davis
% 
% Purpose
% ========
% 
% Script file to test convergence of ground state energy 
% with h-refinement of Spectral FEM basis functions
%
clc;
clear all;
close all

%
% Set timer
%
tstart  = cputime;
% 
global numdofsSch;

% ==========================================
% Problem data (Chemistry part of the input)
% ==========================================
% 
% Set up Problem Parameters:
% 
fprintf('\n');
fprintf('----> Setting up the Problem Parameter . . .\n');
fprintf('\n');

Zion = 6;
qtot = 0.0;             % net charge; nonzero for ionic calculation
nelec = Zion - qtot;

problemdata = struct('Z',Zion, 'nelec',nelec);

% 
% Quadrature
%
global quadtype;
quadtype = 1;          % Type of quadrature: 1 = Gauss, 2 = Lobatto

% 
% 1. Quadrature for radial Schrodinger eq.:
% 
nspSch   = 53;
% nspSch   = pSch + 1;   % Number of quadrature points/weights; 
                       % nsp = p + 1 => diagonal overlap (S)-matrix. This 
                       % leads to exact integration of H-matrix but
                       % "underintegrates" S-matrix. Nonetheless, 
                       % optimal convergence is preserved
    
quaddatainputSch = struct('quadtype',quadtype, 'nsp',nspSch);

% 
% 2. Quadrature for radial Poisson eq.:
% 
nspPois   = nspSch;
% nspPois   = pPois + 1;   % Number of quadrature points/weights; 
                         % nsp = p + 1 => diagonal overlap (S)-matrix. This 
                         % leads to exact integration of H-matrix but
                         % "underintegrates" S-matrix. Nonetheless, 
                         % optimal convergence is preserved
    
quaddatainputPois = struct('quadtype',quadtype, 'nsp',nspPois);

% 
% SCF
% 
maxscfiter  = 200;
scftol      = 1.0E-6;
scftolmax   = 5.0E-5;
mixingparam = 0.7;

scfopts = struct('maxiter',maxscfiter, 'scftol',scftol, 'scftolmax',scftolmax, ...
                 'mixingparam',mixingparam);

maxeigiter = 2000;
eigtol = 1.0E-10;
eigdisp = 0;

eigopts = struct('maxit',maxeigiter, 'tol',eigtol, 'disp',eigdisp);

% 
% Domain:
% 
rmin = 0.0;
rmax = 10.0;

%
% Mesh:
% 
domain      = [rmin rmax];
toplot      = 'on';

% 
% 1. Mesh for radial Schrodinger eq.:
% 
meshtypeSch    = 1;       % type: 1 = uniform, 2 = logarithmic, 3 = exponential
aSch           = 10.0;
numelSch       = 5*(2.^[0:3]);
pSch           = 2;

% 
% 2. Mesh for radial Poisson eq.:
% 
meshtypePois    = meshtypeSch;       % type: 1 = uniform, 2 = logarithmic, 3 = exponential
aPois           = aSch;
pPois           = pSch;
numelPois       = numelSch;

% 
% Initialize:
% 
Etoth  = zeros(length(numelSch),1);
Ekinh  = zeros(length(numelSch),1);
Ecoulh = zeros(length(numelSch),1);
Eenuch = zeros(length(numelSch),1); 
Exch   = zeros(length(numelSch),1);
totdof = zeros(length(numelSch),1); 

for run = 1:length(numelSch)
    meshoptsSch = struct('element_order',pSch, 'numel',numelSch(run), 'domain',domain, ...
                         'verbose',toplot, 'element_type','spectral');
%              
    meshdatainputSch = struct('meshtype',meshtypeSch, 'param',aSch, ...
                              'meshopts',meshoptsSch);

    meshoptsPois = struct('element_order',pPois, 'numel',numelPois(run), 'domain',domain, ...
                          'verbose',toplot, 'element_type','spectral');
%
    meshdatainputPois = struct('meshtype',meshtypePois, 'param',aPois, ...
                               'meshopts',meshoptsPois); 
    
    fprintf('\n');
    fprintf('Solving Radial DFT problem for # of elements = %-5d \n', numelSch(run));
    fprintf('\n')
    
    [lam, u, scf_converged, Etoth(run), Ekinh(run), Ecoulh(run), Eenuch(run), Exch(run), history, exitflag] = radialdft(problemdata, ...
                                                   meshdatainputSch, meshdatainputPois, ...
                                                   quaddatainputSch, quaddatainputPois, ...
                                                   scfopts, eigopts);
    totdof(run) = numdofsSch;
    
    if exitflag ~= 0
        fprintf('\n');
        fprintf('Radial Kohn-Sham Equation could not be solved with Polynomial Order: %-2d and %-5d elements \n', pSch, numelSch(run));
        fprintf('because MATLAB''s eigensolver (eigs) failed! \n');
        fprintf('\n');
        timetot = cputime - tstart;
        fprintf('********************************************************* \n');
        fprintf('\n');
        fprintf(' Total Computing Time so far: %20.3e s \n', timetot);
        fprintf('\n');
        fprintf('********************************************************* \n');
        error('h-Convergence Test Aborted!');    
    else
        fprintf('\n');
        fprintf('Radial Kohn-Sham Equation Successfully solved with Polynomial Order: %-2d and %-5d elements \n', pSch, numelSch(run));
        fprintf('h-Convergence Test Continuing . . . \n');
        fprintf('\n');
    end
    
end

% 
% Plot Total energy:
%     
width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 5.0;     % AxesLineWidth
fsz = 18;      % Fontsize
lw = 2.0;      % LineWidth

figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(totdof, Etoth, 'r', 'LineWidth', 1.25)

% 
xlabel('# of DOF (\equiv # of FE Basis Functions)','Interpreter','tex');
ylabel('Converged Total Energy (Ha)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');

    
timetot = cputime - tstart;
fprintf('*************************************************************** \n');
fprintf('\n');
fprintf(' Time to complete the Simulation = %20.3e s\n', timetot);
fprintf('\n');
fprintf('*************************************************************** \n');

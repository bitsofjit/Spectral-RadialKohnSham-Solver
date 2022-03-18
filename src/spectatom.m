%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Subhajit Banerjee & N. Sukumar
% July 2017
% UC Davis
% 
% Purpose
% ======= 
% 
% Script file to run radialdft:
%cd


clc;
clear all;
close all

addpath(genpath(pwd));
%
% Set timer
%
tstart  = cputime;
% 
% ==========================================
% Problem data (Chemistry part of the input)
% ==========================================
% 
% For benchmarking, check:
% https://www.nist.gov/pml/data/atomic-total-energies-and-eigenvalues-html
% 
% Set up Problem Parameters:
% 
fprintf('\n');
fprintf('----> Setting up the Problem Parameter . . .\n');
fprintf('\n');

Zion = 14;
qtot = 0.0;             % net charge; nonzero for ionic calculation
nelec = Zion - qtot;

problemdata = struct('Z',Zion, 'nelec',nelec);

% 
% Domain:
% 
rmin = 0.0;
rmax = 30.0;

%
% Mesh:
% 
domain      = [rmin rmax];
toplot      = 'on';
% 
% 1. Mesh for radial Schrodinger eq.:
% 
meshtypeSch    = 2;       % type: 1 = uniform, 2 = logarithmic, 3 = exponential
aSch           = 10.0;
pSch           = 6;
numelSch       = 64;
% 
meshoptsSch = struct('element_order',pSch, 'numel',numelSch, 'domain',domain, ...
                     'verbose',toplot, 'element_type','spectral');
%              
meshdatainputSch = struct('meshtype',meshtypeSch, 'param',aSch, ...
                          'meshopts',meshoptsSch);
% 
% 2. Mesh for radial Poisson eq.:
% 
meshtypePois    = 1;       % type: 1 = uniform, 2 = logarithmic, 3 = exponential
aPois           = aSch;
pPois           = 4;
numelPois       = 2*numelSch;
% 
meshoptsPois = struct('element_order',pPois, 'numel',numelPois, 'domain',domain, ...
                      'verbose',toplot, 'element_type','spectral');
%              
meshdatainputPois = struct('meshtype',meshtypePois, 'param',aPois, ...
                           'meshopts',meshoptsPois);                      
% 
% Quadrature
%
global quadtype;
quadtype = 1;          % Type of quadrature: 1 = Gauss, 2 = Lobatto

% 
% 1. Quadrature for radial Schrodinger eq.:
% 
nspSch   = 30;
% nspSch   = pSch + 1;   % Number of quadrature points/weights; 
                       % nsp = p + 1 => diagonal overlap (S)-matrix. This 
                       % leads to exact integration of H-matrix but
                       % "underintegrates" S-matrix. Nonetheless, 
                       % optimal convergence is preserved
%    
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
maxscfiter  = 50;
scftol      = 1.0E-6;
scftolmax   = 5.0E-5;
mixingparam = 0.2;

scfopts = struct('maxiter',maxscfiter, 'scftol',scftol, 'scftolmax',scftolmax, ...
                 'mixingparam',mixingparam);

maxeigiter = 300;
eigtol = 1.0E-10;
eigdisp = 0;

eigopts = struct('maxit',maxeigiter, 'tol',eigtol, 'disp',eigdisp);

[lam, u, scf_converged, Etot, Ekin, Ecoul, Eenuc, Exc, history, exitflag] = radialdft(problemdata, ...
                                                   meshdatainputSch, meshdatainputPois, ...
                                                   quaddatainputSch, quaddatainputPois, ...
                                                   scfopts, eigopts);   
if exitflag ~= 0
    fprintf('\n');
    fprintf('Radial Kohn-Sham Equation could not be solved because \n');
    fprintf('MATLAB''s eigensolver (eigs) failed! \n');
    fprintf('\n');
    timetot = cputime - tstart;
    fprintf('********************************************************* \n');
    fprintf('\n');
    fprintf(' Total Computing Time so far: %20.3e s \n', timetot);
    fprintf('\n');
    fprintf('********************************************************* \n');
    error('spectatom Aborted!');
end

if scf_converged == 1
    fprintf('\n')
    fprintf('       ==========================      \n')
    fprintf('       *** Converged Results  ***      \n')
    fprintf('       ==========================      \n')
    fprintf('\n')
    fprintf('        Etot   = %-11.8e                \n', Etot);
    fprintf('        Ekin   = %-11.8e                \n', Ekin);
    fprintf('        Ecoul  = %-11.8e                \n', Ecoul);
    fprintf('        Eenuc  = %-11.8e                \n', Eenuc);
    fprintf('        Exc    = %-11.8e                \n', Exc);
    fprintf('\n');   
%     fprintf('The ground state energy is: %-11.8e Hartree \n', Etot(end));
else
    fprintf('\n')
    fprintf('SCF did not converged. Maximum number of iterations (%-5d) reached \n', length(history));
    fprintf('\n')
end



    
timetot = cputime - tstart;
fprintf('*************************************************************** \n');
fprintf('\n');
fprintf(' Time to complete the Simulation = %20.3e s\n', timetot);
fprintf('\n');
fprintf('*************************************************************** \n');

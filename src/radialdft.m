function [lam, u, scf_converged, Etot, Ekin, Ecoul, Eenuc, Exc, history, exitflag] = radialdft(problemdata, ...
                                                   meshdatainputSch, meshdatainputPois, ...
                                                   quaddatainputSch, quaddatainputPois, ...
                                                   scfopts, eigopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Subhajit Banerjee & N. Sukumar
% July 2017
% UC Davis
% 
% Purpose
% ======= 
% Solves single-atom, all-electron, radial Kohn-Sham equation by 
% self-consistent field interation (SCF) 
% 
% Inputs
% ====== 
% 
% 
% Outputs
% =======  
% 
% Set up Mesh for Spectral FEM
% 
fprintf('\n');
fprintf('----> Constructing finite element mesh . . . \n');
fprintf('\n');

meshtypeSch = meshdatainputSch.meshtype;
paramSch    = meshdatainputSch.param;  
meshoptsSch = meshdatainputSch.meshopts;
% 
meshsetupSch = mesher(meshtypeSch, paramSch, meshoptsSch);
% 
meshtypePois = meshdatainputPois.meshtype;
paramPois    = meshdatainputPois.param;  
meshoptsPois = meshdatainputPois.meshopts;
% 
meshsetupPois = mesher(meshtypePois, paramPois, meshoptsPois);
% 
% Set up Quadrature
% 
fprintf('\n');
fprintf('----> Setting up quadrature options . . . \n');
fprintf('\n');
% 
quadtypeSch = quaddatainputSch.quadtype;
nspSch      = quaddatainputSch.nsp;
% 
quadsetupSch = getparentquadptswts(quadtypeSch, nspSch);
% 
xiqSch = quadsetupSch.points;
wtqSch = quadsetupSch.weights;
% 
xeSch = meshsetupSch.elemendcoord;
xqSch = getphysicalquadpts(xeSch, xiqSch);
% 
quadtypePois = quaddatainputPois.quadtype;
nspPois      = quaddatainputPois.nsp;
% 
quadsetupPois = getparentquadptswts(quadtypePois, nspPois);
% 
xiqPois = quadsetupPois.points;
wtqPois = quadsetupPois.weights;
% 
xePois = meshsetupPois.elemendcoord;
xqPois = getphysicalquadpts(xePois, xiqPois);

% 
% Initialize density
% 
fprintf('\n');
fprintf('----> Constructing initial density . . . \n');
fprintf('\n');
nelec = problemdata.nelec;
rhoqSch   = initialdensity(nelec, xeSch, xqSch, wtqSch); 
rhoqPois  = initialdensity(nelec, xePois, xqPois, wtqPois);
% 
% Commence SCF Iterations
% 
[lam, u, scf_converged, Etot, Ekin, Ecoul, Eenuc, Exc, history, exitflag] = KSscfiteration(problemdata, rhoqSch, rhoqPois, ...
                                                         meshsetupSch, meshsetupPois, ...
                                                         meshoptsSch, meshoptsPois, ...
                                                         quadsetupSch, quadsetupPois, ...
                                                         xqSch, xqPois, scfopts, eigopts);

if exitflag ~= 0
    fprintf('\n');
    fprintf('SCF Iteration Aborted . . . \n');
    fprintf('\n');
    fprintf('Hit any key to Abort the Radial Kohn-Sham Equation Solver . . . \n');
    pause;
    return
else
    fprintf('\n');
    fprintf('SCF Iterations Ended successfully . . . \n');
    fprintf('\n');
    return
end

end




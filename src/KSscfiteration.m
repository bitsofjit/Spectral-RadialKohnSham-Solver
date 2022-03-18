function [lam, u, scf_converged, Etot, Ekin, Ecoul, Eenuc, Exc, history, eigsflag] = KSscfiteration(problemdata, rhoqSch, rhoqPois,  ...
                                                   meshsetupSch, meshsetupPois, ...
                                                   meshoptsSch, meshoptsPois, ...
                                                   quadsetupSch, quadsetupPois, ... 
                                                   xqSch, xqPois, scfopts, eigopts)                                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Subhajit Banerjee & N. Sukumar
% July 2017
% UC Davis
% 
% Purpose
% ======= 
% Perform SCF iterations:
% Inputs
% ====== 
% 
% 
% Outputs
% =======  

global numdofsSch;

% 
% Unpack:
% 
Zion    = problemdata.Z;
nelec   = problemdata.nelec;
[~, ~, lo, fo] = getatomicstates(Zion);   % no(:), lo(:)  : quantum numbers "n" and "l"
                                          % fo(:)         : occupancy of the (n, l) state
Lmax = max(lo);

xeSch      = meshsetupSch.elemendcoord;
xinSch     = meshsetupSch.parentnodecoord;
connectSch = meshsetupSch.connect;
pSch       = meshoptsSch.element_order;
numelSch   = meshoptsSch.numel;
numnodSch  = pSch*numelSch + 1;
% 
xePois      = meshsetupPois.elemendcoord;
xinPois     = meshsetupPois.parentnodecoord;
connectPois = meshsetupPois.connect;
pPois       = meshoptsPois.element_order;
numelPois   = meshoptsPois.numel;
numnodPois  = pPois*numelPois + 1;
% 
% Set up arrays for DOF related info
% 
fprintf('\n');
fprintf(' ----> Constructing arrays for DOF info \n');
fprintf('\n');
[nodaldofsSch, dofinfoSch, numdofsSch] = KSsetupdofs(numnodSch);
[nodaldofsPois, dofinfoPois, numdofsPois] = Poissetupdofs(numnodPois);

xiqSch = quadsetupSch.points;
wtqSch = quadsetupSch.weights;

xiqPois = quadsetupPois.points;
wtqPois = quadsetupPois.weights;
% 
% Set up basis functions in parent domain
% 
fprintf('\n');
fprintf('----> Constructing basis functions . . . \n');
fprintf('\n');
[phiqSch, dphiqSch] = getshapefunc(pSch, xiqSch);
[phiqPois, dphiqPois] = getshapefunc(pPois, xiqPois);

maxiter     = scfopts.maxiter;
scftol      = scfopts.scftol;
scftolmax   = scfopts.scftolmax;
mixingparam = scfopts.mixingparam;


% Initialization
% lam = zeros(neigs,Lmax+1);
% u = zeros(numdofs,neigs,Lmax+1);  

lam = cell(Lmax+1,1);        % to store eigenvalues of H u = lam S u:
                             % lam(i,j) = i-th eigenvalue of j-th 
                             % L-value, j=0..Lmax
u = cell(Lmax+1,1);          % to store eigenvectors of H u = lam S u;
                             % u(i,j,k) = i-th component of j-th 
                             % vector of k-th L-value, k=0..Lmax

% Computing the number of eigenvalues required as each l = 0, ..., Lmax
neigs = zeros(Lmax + 1,1);
for l = 0:Lmax
    neigs(l+1) = length(find(lo == l));
end

scf_converged = false;

% SCF iterations begins
fprintf('\n');
fprintf('----> SCF iterations begin . . . \n');
fprintf('\n');

for iter = 1:maxiter       
    fprintf('\n');
    fprintf(' ----> Starting SCF Iteration # %-5d . . .\n', iter);
    fprintf('\n');
    %     
    % Normalize density on Poisson and K-S mesh:
    %     
    rhointSch  = sphericalintegral(xeSch, xqSch, wtqSch, rhoqSch);
    rhointPois = sphericalintegral(xePois, xqPois, wtqPois, rhoqPois);
    
    fprintf('\n');
    fprintf(' Integral of input density before normalization (on K-S Mesh): %-11.8e \n', rhointSch);
    fprintf('\n');
    fprintf(' Integral of input density before normalization (on Poisson Mesh): %-11.8e \n', rhointPois);
    fprintf('\n');
    
    rhoqSch = -rhoqSch./rhointSch*nelec;        % charge-density on K-S Mesh
    rhoqPois = -rhoqPois./rhointPois*nelec;     % charge-density on Poisson Mesh
    
    % Construct Potential:
    fprintf('\n');
    fprintf('----> Constructing Hartree potential . . . \n');
    fprintf('\n');
    % First get nodal solution for Vp(r) = r * Vh(r)
    % -ive sign as Vp" = -4 * pi * n(r) * r
    VpPois = solveradialpoisson(xePois, connectPois, numelPois, ...
                                nodaldofsPois, dofinfoPois, numdofsPois, ...
                                xqPois, wtqPois, phiqPois, dphiqPois, -4.0*pi*rhoqPois);
                               
    % Transform to K-S quadrature grid 
    VpqSch = fetoquadgen(xinPois, xiqSch, xePois, xeSch, connectPois, ...
                         [0.0; VpPois; 0.0]);               % zero Dirichlet on the left 
                                                            % and zero Dirichlet on the right      
    % Construct solution satisfying the BC Vh(r=rmax) = Z/rmax, and we 
    % divide by "r" to obtain the Hartree potential.
    % This is because we were solving for Vp(r) = r * Vh(r):
    VhqSch = VpqSch./xqSch + Zion/xeSch(end);
    
    
    % Vnq(:,:)      : quadrature-grid  representatin of Vn(x) = -Z/(|x - R|);
    %                 Vn(i,j) = value at i-th quadrature point of j-th element
    %Vnq = computenuclearpotential(Zion, xqSch); 
    
    %Veffq = Vhq + Vnq;
    
    VeffqSch = VhqSch - Zion./xqSch;        % Add Nuclear potential
    
    [excqSch, VxcqSch] = exhangecorrelation('lda', -rhoqSch); 
    
    VeffqSch = VeffqSch + VxcqSch;
    
    % Find eigenstates
    fprintf('\n');
    fprintf('----> Computing eigenstates . . . \n');
    fprintf('\n');
    
    for l = 0:Lmax
        % u{l+1,1} = zeros(numdofs, neigs(l+1));
        % lam{l+1,1} = zeros(neigs(l+1),1);
        [lam{l+1,1}, u{l+1,1}, eigsflag] = solveradialschroed(xeSch, connectSch, numelSch, ...
                                           nodaldofsSch, dofinfoSch, numdofsSch, ...
                                           VeffqSch, l, xqSch, xiqSch, wtqSch, phiqSch, ...
                                           dphiqSch, neigs(l+1), eigopts); 
        if eigsflag ~= 0
            break
        else 
            continue
        end
    end
    
    if eigsflag ~= 0
        scf_converged = false;
        Etot  = NaN; 
        Ekin  = NaN;
        Ecoul  = NaN;
        Eenuc  = NaN; 
        Exc  = NaN; 
        history  = NaN;
        return
    end
    
    % Construct density:
    fprintf('\n');
    fprintf('----> Constructing output density . . . \n');
    fprintf('\n');
    
    rhooutqSch = 0.0;
    rhooutqPois = 0.0;
    Eband = 0.0;
    
    for l = 0:Lmax
        nbands = neigs(l+1);
        U = reshape([u{l+1}(:,:)], [numdofsSch,nbands]);
        LAM = reshape([lam{l+1,1}(:)], [nbands, 1]);
        occ = fo(lo == l);
        for n = 1:nbands
            % transform to quadrature grid
            fullc = zeros(numdofsSch+2,1);
            fullc(2:numdofsSch+1) = U(:,n);
            
            % Generate quadrature point values of P_{nl} on Schrodinger mesh:
            PqSch = fetoquad(xinSch, xiqSch, connectSch, fullc);
            
            % Generate quadrature point values of P_{nl} on Poisson mesh:
            PqPois = fetoquadgen(xinSch, xiqPois, xeSch, xePois, connectSch, fullc);
            %
            % Form Charge density and accumulate
            %
            % On Schrodinger mesh:
            rhooutqSch = rhooutqSch - occ(n)*PqSch.^2./(xqSch.^2)/(4.0*pi);
            
            % On Poisson mesh:
            rhooutqPois = rhooutqPois - occ(n)*PqPois.^2./(xqPois.^2)/(4.0*pi);
            
            Eband = Eband + occ(n)*LAM(n);
            
        end
    end
 
    
    [T_s(iter), E_ee(iter), E_en(iter), E_xc(iter), E_tot(iter)] = ...
                computetotalenergy(xeSch, xqSch, wtqSch, rhooutqSch, VeffqSch, VhqSch, ...
                                   excqSch, Zion, Eband);
    history(iter) = iter;
    
    rhointSch = sphericalintegral(xeSch, xqSch, wtqSch, rhooutqSch);
    rhointPois = sphericalintegral(xePois, xqPois, wtqPois, rhooutqPois);
    
    fprintf('\n');
    fprintf(' Integral of OUTPUT density before normalization: %-11.8e \n', rhointSch);
    fprintf('\n');
    
    rhooutqSch  = -rhooutqSch./rhointSch*nelec;
    rhooutqPois = -rhooutqPois./rhointPois*nelec;
   
    % Assess convergence, exit if achieved
    fprintf('\n');
    fprintf('----> Assessing convergence . . . \n');
    fprintf('\n');
    
    drhoq = norm((rhooutqSch - rhoqSch),'fro')/norm(rhoqSch,'fro');
    
    fprintf('\n');
    fprintf('Fractional difference of input and output densities = %-11.8e \n', drhoq);
    fprintf('\n');

    if (iter > 3)
        Emax = max(E_tot(iter-3:iter));
        
        Emin = min(E_tot(iter-3:iter));
        
        if (drhoq < scftol && Emax-Emin < scftolmax)
            fprintf('\n');
            fprintf('Self-consistent convergence criterion achieved \n');
            fprintf('\n');
            fprintf('SCF Iteration: %-4d \n', iter);
            fprintf('tolscf: %-11.8e \n', scftol);
            fprintf('drhoq: %-11.8e \n', drhoq);
            fprintf('tolscfmax: %-11.8e \n', scftolmax);
            fprintf('Emax - Emin: %-11.8e \n', Emax-Emin);
            fprintf('\n');
            
            scf_converged = true;
            break
        end
    end
    
    fprintf('\n');
    fprintf('----> SCF iterations end . . . \n');
    fprintf('\n');
          
   % if not converged, form new input density
   scf_converged = false;
   fprintf('\n');
   fprintf('----> Forming new input density on Poisson mesh. . . \n');
   fprintf('\n');
   
   rhoqPois =  linearmix(rhoqPois, rhooutqPois, mixingparam);
   rhoqSch  =  linearmix(rhoqSch, rhooutqSch, mixingparam);

end

Etot  = E_tot(end);
Ekin  = T_s(end);
Ecoul = E_ee(end);
Eenuc = E_en(end);
Exc   = E_xc(end);
    
end
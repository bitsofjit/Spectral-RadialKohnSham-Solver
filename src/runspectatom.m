function runspectatom(problemdata, meshdatainputSch, meshdatainputPois, ...
                      quaddatainputSch, quaddatainputPois, scfopts, eigopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Subhajit Banerjee & N. Sukumar
% July 2017
% UC Davis
% 
% Purpose
% ======= 
% 
%

addpath(genpath(pwd));
        

[lam, u, scf_converged, Etot, Ekin, Ecoul, Eenuc, Exc, history] = radialdft(problemdata, ...
                                                   meshdatainputSch, meshdatainputPois, ...
                                                   quaddatainputSch, quaddatainputPois, ...
                                                   scfopts, eigopts);   

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

% % 
% % Plot Total energy:
% %     
% width = 6;     % Width in inches
% height = 4;    % Height in inches
% alw = 5.0;     % AxesLineWidth
% fsz = 18;      % Fontsize
% lw = 2.0;      % LineWidth
% 
% figure
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
% set(gca, 'FontSize', fsz, 'LineWidth', alw);
% set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
% hold on;
% % 
% plot(history, Etot, 'r', 'LineWidth', 1.25)
% 
% % 
% xlabel('Iteration Number','Interpreter','tex');
% ylabel('Total Energy (Ha)','Interpreter','tex');
% box on; 
% axis square;
% set(gca,'LineWidth',lw);
% set(gca,'FontSize',fsz);
% set(gca,'FontWeight','Bold');
% set(gcf,'color','w');

    
timetot = cputime - tstart;
fprintf('*************************************************************** \n');
fprintf('\n');
fprintf(' Time to complete the Simulation = %20.3e s\n', timetot);
fprintf('\n');
fprintf('*************************************************************** \n');

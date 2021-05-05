%% Steady State: Compare results of Chebfun, (Semi-) Anayltical and OpenFOAM solutions (Infinite/ Periodic Fluid Channel)
function [foamDataYCPV,foamDataYU,chebC,chebP,chebU,chebV] = infChannelCompare(foamEndTime,plots)
    close all;
%% Initialize shared variables
    % Pressure figure
    a = subplot(2,2,1);%figure(1);
    title("Pressure $p$",'Interpreter','latex');
    xlabel("y");
    % Velocity figure
    b = subplot(2,2,2);%figure(2);
    title("Fluid velocity $u_1$",'Interpreter','latex');
    xlabel("y");
    % Concentration figure
    c = subplot(2,2,3);%figure(3);
    title("Concentration $c$",'Interpreter','latex');
    xlabel("y");
    % Electric Potential figure
    d = subplot(2,2,4);%figure(4);
    title("Electric Potential $\phi$",'Interpreter','latex');
    xlabel("y");
%% Import OpenFOAM post-process files/data
    
    % File names (constant for all times)
    foamFileNames = ["line_C.anions_p_V.xy" "line_U.xy"];
    
    % Data storage (End/Steady time) for y-pos, concentration, pressure and
    % potential
    foamDataYCPV = load(strcat("../tutorials/pnpFoam/planePoiseuille_cyclic/postProcessing/singleGraph/", ...
        num2str(foamEndTime),"/",foamFileNames(1)));
    % Data storage (End/Steady time) for y-pos, velocity
    foamDataYU = load(strcat("../tutorials/pnpFoam/planePoiseuille_cyclic/postProcessing/singleGraph/", ...
        num2str(foamEndTime),"/",foamFileNames(2)));
    % y-pos storage (same in both openFOAM files due to sampling type the same)
    foamYPos = foamDataYU(:,1);
%% Plot OpenFOAM data

    if plots
        % Plot pressure
        subplot(a); hold on;%figure(a); hold on;
        plot(foamYPos, foamDataYCPV(:,3),'r'); hold off;
        % Plot velocity
        subplot(b); hold on;%figure(b); hold on;
        plot(foamYPos,foamDataYU(:,2),'r'); hold off;
        % Plot concentration
        subplot(c); hold on;%figure(c); hold on;
        plot(foamYPos, foamDataYCPV(:,2),'r'); hold off;
        % Plot electric potential
        subplot(d); hold on;%figure(d); hold on;
        plot(foamYPos, foamDataYCPV(:,4),'r'); hold off;
    end
%% Run Chebfun scripts and import results

    % !!! Currently have to change physical parameter vals in script manually !!!
    [chebC,chebP,chebU,chebV] = infChannelCheb(false,false);
    
%% Plot Chebfun data

    if plots
        % Plot pressure
        subplot(a); hold on;%figure(a); hold on;
        plot(chebP,'b*'); legend('OpenFOAM','Chebfun','Interpreter','latex'); hold off;
        % Plot velocity
        subplot(b); hold on;%figure(a); hold on;
        plot(chebU,'b*'); legend('OpenFOAM','Chebfun','Interpreter','latex'); hold off;
        % Plot concentration
        subplot(c); hold on;%figure(a); hold on;
        plot(chebC,'b*'); legend('OpenFOAM','Chebfun','Location','NorthWest','Interpreter','latex'); hold off;
        % Plot electric potential
        subplot(d); hold on;%figure(a); hold on;
        plot(chebV,'b*'); legend('OpenFOAM','Chebfun','Location','NorthWest','Interpreter','latex'); hold off;
    end
    
% %% Run (Semi-) Analytical/ Asymptotic scripts and import results
% 
%     % !!! Currently have to change physical parameter vals in script manually %
%     [asymY,asymC,asymV,~,~] = infChannelAnaSol(false,false);
%     
% %% Plot (Semi-) Analytical/ Asymptotic data
%     % Plot pressure
%     %plot(asymY,asymP,'b'); !!! TO DO !!!!!! TO DO !!!!!! TO DO !!!!!! TO DO !!!
%     % Plot velocity
%     %plot(asymY,asymU,'b'); !!! TO DO !!!!!! TO DO !!!!!! TO DO !!!!!! TO DO !!!
%     % Plot concentration
%     figure(c); hold on;
%     plot(asymY,asymC,'r'); hold off;
%     % Plot electric potential
%     figure(d); hold on;
%     plot(asymY,asymV,'r'); hold off;
%     
%     
%     
%     
%     
    
    
    
    

end
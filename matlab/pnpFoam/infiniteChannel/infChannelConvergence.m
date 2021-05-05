%% Construct L2 norm between Chebfun and OpenFOAM for infinite Channel over increasing refinement levels
function []=infChannelConvergence(foamEndTime)
close all
%% Run infChannelCheb for Chebfun data

% !!! Currently have to change physical parameter vals in script manually !!!
    [chebC,chebP,chebU,chebV] = infChannelCheb(false,false);
    
%% Variables
 
    % File names (constant for all times)
    foamFileNames = ["line_C.anions_p_V.xy" "line_U.xy"];
    
    % refinement levels (# of cells in y)
    numCells = [20,40,80,160,320];
    
    % L2 norm storage at refinement levels (row:= Field, col:= ref lvl)
    % row 1 = Concentration L2norm, row 2 = Pressure L2norm, row 3 =
    % Potential L2norm
    foamL2Norm = zeros(3,length(numCells));
    
%% Load data and compute L2 norm
    % Loop through # of refinements
    for j =1:length(numCells)

        % Pause until the data is available (press any key while focused on MatLab when ready)
        disp("Run openFOAM for data files, press any key when ready");
        pause;
        
        % Load y-pos, concentration, pressure and potential
        foamDataYCPV = load(strcat("../tutorials/pnpFoam/planePoiseuille_cyclic/postProcessing/singleGraph/", ...
            num2str(foamEndTime),"/",foamFileNames(1)));
        
        % Calculate L2norm for concentration
        foamL2Norm(1,j) = norm(foamDataYCPV(:,2) - chebC(foamDataYCPV(:,1)),2)./norm(chebC(foamDataYCPV(:,1)),2);
        
        % Calculate L2norm for pressure
        foamL2Norm(2,j) = norm(foamDataYCPV(:,3) - chebP(foamDataYCPV(:,1)),2)./norm(chebP(foamDataYCPV(:,1)),2);
        
        % Calculate L2norm for potential
        foamL2Norm(3,j) = norm(foamDataYCPV(:,4) - chebV(foamDataYCPV(:,1)),2)./norm(chebV(foamDataYCPV(:,1)),2);
    end
    
%% Plot L2norm against refinement lvls
figure;
loglog(numCells,foamL2Norm(1,:),'r*-'); hold on;
loglog(numCells,foamL2Norm(2,:),'b*-');
loglog(numCells,foamL2Norm(3,:),'k*-');
loglog(numCells,numCells.^-1,'g');
legend('c','$\phi$','p','$\mathcal{O}(Number Cells ^{-1})$','interpreter','latex');
title("L2norm convergency plot: Infinite Channel");
xlabel("Refinement level (number of cells in y)");
ylabel("L2norm");

end
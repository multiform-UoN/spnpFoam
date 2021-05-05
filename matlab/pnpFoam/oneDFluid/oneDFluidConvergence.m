%% Provides a convergency plot of pnpFoam, taking the L2norm of the difference between Chebfun and OpenFOAM

%% Initialize variables

% File name (constant for all times)
file_name = "/line_C.anions_C.cations_V.xy";

% End time
t_end = 1e-05;

% vector of # of cells
numCells = [100,1000,5000,10000];

% vectors of L2 norms
L2normsCaCcV = zeros(3,length(numCells));

%% Compute 'exact' or 'true' result (pnp.m tol=1e-12)
[t,chebCc,chebCa,chebV] = oneDFluidCheb(1e-11,false);

%% Load data and compute L2 norm
    % Loop through # of refinements
    for j =1:length(numCells)
        
        % Pause until the data is available (press any key while focused on MatLab when ready)
        disp("Run openFOAM for data files, press any key when ready");
        pause; 
        
            % Load data
            foamDataXCaCcV = load(strcat("../../../tutorials/pnpFoam/testCase1D/postProcessing/singleGraph/",num2str(t_end),file_name));
            
            % Compute L2 norm (compute cheb poly's at FOAM points)
            L2normsCaCcV(1,j) = norm(foamDataXCaCcV(:,2) - chebCa(foamDataXCaCcV(:,1),end),2)./norm(chebCa(foamDataXCaCcV(:,1),end),2);
            L2normsCaCcV(2,j) = norm(foamDataXCaCcV(:,3) - chebCc(foamDataXCaCcV(:,1),end),2)./norm(chebCc(foamDataXCaCcV(:,1),end),2);
            L2normsCaCcV(3,j) = norm(foamDataXCaCcV(:,4) - chebV(foamDataXCaCcV(:,1),end),2)./norm(chebV(foamDataXCaCcV(:,1),end),2);
    end

%% Plot L2 norm against # of cells
figure;
loglog(numCells,L2normsCaCcV(1,:),'r*-','LineWidth',2);
hold on;
loglog(numCells,L2normsCaCcV(2,:),'b*-','LineWidth',2);
loglog(numCells,L2normsCaCcV(3,:),'k*-','LineWidth',2);
loglog(numCells,numCells.^-1,'g','LineWidth',2);
legend('$c_{2}$','$c_{1}$','$\phi$','interpreter','latex');
xlabel('N (number of cells in x)','interpreter','latex');
title('$L^{2}$ Convergency plot: Single Region','interpreter','latex');
set(gca,'FontSize',24);


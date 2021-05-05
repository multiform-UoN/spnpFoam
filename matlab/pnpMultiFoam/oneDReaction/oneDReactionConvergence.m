%% Provides a convergency plot for multiRegion, taking the L2norm of the difference between Chebfun and OpenFOAM

%% Initialize variables

% File name (constant for all times)
fileName = "/line_T.xy";
% Directory for fluid side
fluidFileDir = "singleGraphFluid/fluid/";
% Directory for solid side
solidFileDir = "singleGraphSolid/solid/";

% End time(s)
t_end = [2,2,3,3,6,18,57,208];

% vector of # of cells (3:2 split between fluid:solid)
numCells = [10,20,40,80,160,320,640,1280];

% vectors of L2 norms of T for varying refinement levels
L2normsT = zeros(1,length(numCells));
%% Compute 'exact' or 'true' result, i.e. Chebfun results
[chebTFluid, chebTSolid] = oneDReactionChebSteady(1e-12,false);

%% Load data and compute L2 norm
    % Loop through # of refinements
    for j =1:length(numCells)

        % Pause until the data is available (press any key while focused on MatLab when ready)
        disp("Run openFOAM for data files, press any key when ready");
        pause; 

            % Load fluid side data
            line_T_Fluid = load(strcat("~/OpenFOAM/multiformFoam/tutorials/transport/scalarMultiRegionFoam/simpleOneD/postProcessing/", ...
                fluidFileDir,num2str(t_end(j)),fileName));
            % Load solid side data
            line_T_Solid = load(strcat("~/OpenFOAM/multiformFoam/tutorials/transport/scalarMultiRegionFoam/simpleOneD/postProcessing/", ...
                solidFileDir,num2str(t_end(j)),fileName));
            % Join solid and fluid data together
            line_T_Both = [line_T_Fluid ; line_T_Solid];
            % Compute L2 norm of anions (compute cheb poly's at FOAM points)
            L2normsT(1,j) = norm(line_T_Both(:,2) - [chebTFluid(line_T_Fluid(:,1),end) ; chebTSolid(line_T_Solid(:,1),end)],2)./norm([chebTFluid(line_T_Fluid(:,1),end) ; chebTSolid(line_T_Solid(:,1),end)],2);

    end

%% Plot L2 norm against # of cells
figure;
loglog(numCells,L2normsT(1,:),'r*-','LineWidth',2);
hold on;
loglog(numCells,numCells.^(-1),'g','LineWidth',2)
xlabel('N (number of cells in x)','interpreter','latex');
legend('$c_{i}$','$\mathcal{O}(N^{-1})$','interpreter','latex');
set(gca,'fontSize',24);


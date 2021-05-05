%% Plots data from OpenFOAM & Chebfun into Matlab figures (multiRegion, reactive mappedChemicalKinetics)

%% Initialize variables
% File name (constant for all times)
fileName = "/line_T.xy";
% Directory for fluid side
fluidFileDir = "singleGraphFluid/fluid/";
% Directory for solid side
solidFileDir = "singleGraphSolid/solid/";
% Initial time
iniTime = 0.5;
% Time step
dt = 0.5;
% End time
endTime = 16.5;
% # of time steps
steps = endTime/dt;
% Times plotted + end time
interSteps = [iniTime,2,5,endTime];
% Create figure
figure;
hold on;
%% Load OpenFOAM data and plot
% Loop through all time directories in interSteps and plot
for i=1:length(interSteps)
    
    % Load fluid side data
    line_T_Fluid = load(strcat("~/OpenFOAM/multiformFoam/tutorials/scalarMultiRegionFoam/simpleOneD/postProcessing/", ...
        fluidFileDir,num2str(interSteps(i)),fileName));
    % Load solid side data
    line_T_Solid = load(strcat("~/OpenFOAM/multiformFoam/tutorials/scalarMultiRegionFoam/simpleOneD/postProcessing/", ...
        solidFileDir,num2str(interSteps(i)),fileName));
    
    % Plot T fluid
    if i == length(interSteps)
        p(i+1)= plot(line_T_Fluid(:,1),line_T_Fluid(:,2),'r--','LineWidth',2);
    else
        p(i+1)= plot(line_T_Fluid(:,1),line_T_Fluid(:,2),'k--','LineWidth',2);
    end

    % Plot T solid
    hold on;
    if i== length(interSteps)
        k(i+1)=plot(line_T_Solid(:,1),line_T_Solid(:,2),'r','LineWidth',2);
    else
        k(i+1)=plot(line_T_Solid(:,1),line_T_Solid(:,2),'k','LineWidth',2);
    end

end

%% Plot t=0 lines (Can't obtain from openFOAM?!)

% Plot t=0 for fluid side
p(1) = plot(line_T_Fluid(:,1),exp(-200.*(line_T_Fluid(:,1)+1).^2),'b--','LineWidth',2);

% Plot t=0 for solid side
k(1) = plot(line_T_Solid(:,1),zeros(length(line_T_Solid(:,1)),1),'b','LineWidth',2);

%% Plot transient case of Chebfun (linear times)
%Load Chebfun data
[ChebTFluid, ChebTSolid] = oneDReactionChebTransient(1e-11,false,true);
% Intermediate time step Chebfun indexes
ChebInterTimes = interSteps(1:end-1)./dt + 1;
% Chebfun initial time index
ChebInitTime = 1;
% Chebfun final time index
ChebFinalTime = endTime./dt +1;
% Plot values at specified time indexes
chInter = plot(ChebTFluid(:,ChebInterTimes),'k*',ChebTSolid(:,ChebInterTimes),'k*');
chInit = plot(ChebTFluid(:,ChebInitTime),'b*',ChebTSolid(:,ChebInitTime),'b*');
chFinal = plot(ChebTFluid(:,ChebFinalTime),'r*',ChebTSolid(:,ChebFinalTime),'r*');

%% Style the plots
xlabel('x [m]','interpreter','latex');
ylabel('Concentration [$\frac{mol}{m}$]','interpreter','latex');
legend([p(1) k(1) p(2) k(2) p(end) k(end)],'$c_{f}(0)$','$c_{s}(0)$', ...
    ['$c_{f}(0<t<$' num2str(endTime) ')'],['$c_{s}(0<t<$' num2str(endTime) ')'],['$c_{f}($' num2str(endTime) ')'],['$c_{s}($' num2str(endTime) ')'],'interpreter','latex');
set(gca,'FontSize',30);

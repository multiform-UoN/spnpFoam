%% Post-processing of twoDRandomField case for pnpMultiFoam
% Takes n simulation runs of pnpMultiFoam across generated random
% porous solid-fluid meshes

% Partitioning of simulations:
% For P porosities take S random seeds and R refinement lvls
% Therefore, n = P*S*R
% Example: For n=27, take P=3, S=3, R=3
%% Variables
% Number of porosities
P = 3;
% Number of refinement lvls
R = 3;
% Number of random seeds
S = 3;

% Array of data files (size SxRxP, i.e an SxR matrix for every P 'page')
fileNameArray = strings(S,R,P);

% Charge numbers
zNa = 1 ; zCl = -1;
%% Populate fileDataArray per porosity 'page'
    for i=1:P
        PoroString = strcat("Poro",num2str(i));
        for j=1:S
            SeedString = strcat("Seed",num2str(j));
            for k=1:R
                RefineString = strcat("Refine",num2str(k),".dat");
                fileNameArray(j,k,i) = strcat(PoroString,SeedString,RefineString);
            end
        end
    end
%% Create figures
% Volume integrals of concentrations
vol = figure(1);
% Total charge density
chrg = figure(2);
% Line color (for porosity)
line_Color = ['b','g','r'];
% Line style (for seed)
line_Style = ['-','--',':'];
% Mark style (for refinement)
line_Marker = ['*','.','o'];
%% Extract data and add to plots
% Per fileDataArray element, extract data (volume integrals + calculate tot charge density)

% Not sure how to do the legend as that would require 27 items. Possibly discuss it in the caption
% instead?
    for i=1:P
        for j=1:S
            for k=1:R
                Style = strcat(line_Color(i),line_Style(j),line_Marker(k));
                % Extract data from relevant file
                fileDataVolInts_tCClVolCNa = importdata(fileNameArray(j,k,i),' ',4);
                figure(vol);
                hold on;
                % Plot tot conc of C.Cl
                plot(fileDataVolInts_tCClVolCNa(:,1),fileDataVolInts_tCClVolCNa(:,2),Style);
                % Plot tot conc of C.Na
                plot(fileDataVolInts_tCClVolCNa(:,1),fileDataVolInts_tCClVolCNa(:,3),Style);
                hold off;
                figure(chrg);
                hold on;
                % Compute tot chrg density = tot(C.Na)*zNa + tot(C.Cl)*zCl
                fileDataVolInt_ChrgDensity = [fileDataVolInts_tCClVolCNa(:,1), ...
                    fileDataVolInts_tCClVolCNa(:,2).*zCl + fileDataVolInts_tCClVolCNa(:,3).*zNa];
                % Plot tot chrg density
                plot(fileDataVolInt_ChrgDensity(:,1),fileDataVolInt_ChrgDensity(:,2),Style);
                hold off;
            end
        end
    end


function oneDFluidFoam(intTimes)

    %% Plots data coming from OpenFOAM (singleRegion)

    close all;
    
    %% Initialize variables
    % Time directory folder name (initialize to first time)
    file_dir = 1e-06;
    % File name (constant for all times)
    file_name = "/line_C.anions_C.cations_V.xy";
    % Time step
    dt = 1e-06;
    % End time
    t_end = 1e-05;
    % # of time steps to load
    steps = round(t_end/dt);
    % Create figures
    figure(1);
    hold on;
    figure(2);
    hold on;
    %% Compute Chebfun data
    [t, c1, c2, v] = oneDFluidCheb(1e-11,false,false);
    %% Load data and plot
    % Loop through all time directories
    for i=1:steps

        % Load data
        line_C_anions_C_cations_V = load(strcat("../tutorials/pnpFoam/testCase1D/postProcessing/singleGraph/",num2str(file_dir),file_name));
        % Change figure
        figure(1);
        % Plot Anions

        %subplot(2,1,2);
        %hold on;
        if i == steps
            p(i+1)= plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,2),'r-','LineWidth',3);
        else
            if intTimes==true
                p(i+1)= plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,2),'b-','LineWidth',2);%,'MarkerSize',10);
            end
        end
        hold on;
        % Plot Cations
        %subplot(2,1,1);
        %hold on;
        if i== steps
            k(i+1)=plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,3),'r--','LineWidth',3);
        else
            if intTimes==true
                k(i+1)=plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,3),'b--','LineWidth',2);%,'MarkerSize',10);
            end
        end
        % Change figure
        figure(2);
        hold on;
        % Plot Potential
        if i== steps
            r(i+1)=plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,4),'r-','LineWidth',3);
        else
            if intTimes==true
                r(i+1)=plot(line_C_anions_C_cations_V(:,1),line_C_anions_C_cations_V(:,4),'b--','LineWidth',2);%,'MarkerSize',20);
            end
        end

        % Update time
        file_dir = dt + file_dir;
    end

    %% Plot t=0 lines (Can't obtain from openFOAM?!)
    figure(1);

    hold on;
    p(1) = plot(line_C_anions_C_cations_V(:,1),ones(length(line_C_anions_C_cations_V(:,1)),1).*1e-03,'k-','LineWidth',2);
    l(1) = plot(c1(:,end),'b*','MarkerSize',10);
    l(2) = plot(c2(:,end),'b.','MarkerSize',10);
    k(1) = plot(line_C_anions_C_cations_V(:,1),ones(length(line_C_anions_C_cations_V(:,1)),1).*1e-03,'k--','LineWidth',2);

    figure(2);
    hold on;
    r(1) = plot(line_C_anions_C_cations_V(:,1),0.1.*(1-cos(pi.*line_C_anions_C_cations_V(:,1)./1e-06))./2,'k-','LineWidth',2);
    m(1) = plot(v(:,end),'b*');
    %% Style the plots
    figure(1);

    %subplot(2,1,2);

    title("OpenFOAM: Molar Concentrations $c_{i}$",'interpreter','latex');
    if intTimes==true
        legend([p(2) p(end) p(1) k(2) k(end) k(1)],['0 < t < ' num2str(t_end)],['t = ' num2str(t_end)],'t = 0'...
            ,['0 < t < ' num2str(t_end)],['t = ' num2str(t_end)],'t = 0','Location','northwest','interpreter','latex');
    else
        legend([p(end) p(1) k(end) k(1)],['$c_{2}($' num2str(t_end) '$)$'],'$c_{2}(0)$'...
            ,['$c_{1}($' num2str(t_end) '$)$'],'$c_{1}(0)$','Location','northwest','interpreter','latex');
    end
    xlabel('x [m]','interpreter','latex');
    ylabel(['Concentration ' '[$\frac{mol}{m}$]'],'interpreter','latex');
    set(gca,'FontSize',24);

    %subplot(2,1,1);

    %title("OpenFOAM: Molar Concentration $c_{cations}$",'interpreter','latex');
    %legend([k(2) k(end) k(1)],['0 < t < ' num2str(t_end)],['t = ' num2str(t_end)],'t = 0');
    %xlabel('Domain (x)','interpreter','latex');
    %set(gca,'FontSize',18);

    figure(2);
    %hold on;
    title("OpenFOAM: Electric Potential $\phi$",'interpreter','latex');
    if intTimes==true
        legend([r(2) r(end) r(1)],['0 < t < ' num2str(t_end)],['t = ' num2str(t_end)],'t = 0','Location','southeast');
    else
        legend([r(end) r(1)],['t = ' num2str(t_end)],'t = 0','Location','southeast','interpreter','latex');
    end
    xlabel('x [m]','interpreter','latex');
    ylabel(['Potential ' '[V]'],'interpreter','latex');
    ylim([0,0.1]);
    set(gca,'FontSize',24);
end
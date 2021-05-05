function [t, c1, c2, v] = oneDFluidCheb(tol,plotGraphs,intTimes)

    %% pnp.m -- Poisson Nernst-Planck in 1D

    
    %close all;

    %% Time and space domain
    % Create an interval of the space domain...
    dom = [0 1];
    t_end = 10^(-5);% 10^-5
    dt = t_end/100;
    %...and specify a sampling of the time domain:
    t = 0:dt:t_end;

    %% Physical parameters
    c0 = 1e-3; % ref concentration
    phi0 = 0.1; % ref potential
    L0 = 1e-6; % ref length
    vel=0.000;   % advection flux
    vel0 = 0; % Velocity at walls
    D = 1e-6;   % diffusion
    epsilon = 80*8.8541878176e-12; % electrical permeability
    ele_charge = 1.6021766208e-19; % elementary charge
    boltz_const = 1.38064852e-23; % Boltzmann constant
    Ava_const = 6.022140857e+23; % Avogadro constant
    F = ele_charge*Ava_const; % Faraday constant
    T = 300; % Temperature, set constant
    Z1 = 1; %Charge type
    Z2 = -1; %Charge type
    npc1 = (D*Z1*ele_charge)/(boltz_const*T); % NP flux coeff
    npc2 = (D*Z2*ele_charge)/(boltz_const*T); % NP flux coeff
    % kappa = eta*10; % potential-induced flux
    % potential and concentration are normalised between 0 and 1

    %% PDEs
    % Make the right-hand side of the PDE.
    dom = dom*L0;
    pdefun = @(t,x,c1,c2,v) [...
            diff(D.*diff(c1)+c1.*(npc1.*diff(v)-vel)); ... % Nernst-Planck (c1), (replace minus for vel component)
            diff(D.*diff(c2)+c2.*(npc2.*diff(v)-vel)); ... % Nernst-Planck (c2),
            diff(epsilon.*diff(v))+F.*Z1.*c1 + F.*Z2.*c2...  % Poisson
            ];
    pdeflag = [1  1  0]; % Zero when a variable is indep of time.

    % Assign boundary conditions.
    % No flux and fixed potential
    bc.left  = @(t,c1,c2,v) [D*diff(c1)+c1.*(npc1*diff(v)-vel0); ... % No NP flux (c1)
                          D*diff(c2)+c2.*(npc2*diff(v)-vel0) ; v]; % No NP flux (c2) + fixed potential
    bc.right = @(t,c1,c2,v) [D*diff(c1)+c1.*(npc1*diff(v)-vel0); ...  % No NP flux (c1)
                          D*diff(c2)+c2.*(npc2*diff(v)-vel0); v-phi0]; % No NP flux (c2) + fixed potential

    % Construct a chebfun of the space variable on the domain,
    x = chebfun(@(x) x, dom);
    % and of the initial conditions.
    c1_0 = c0*(x*0+1);
    %c1_0 = c0*(-cos(pi*x/dom(2))+1)/2;
    c2_0 = c0*(x*0 +1);
    %c2_0 = c0*(-cos(pi*-x/dom(2))+1)/2;
    v0 = phi0*(-cos(pi*x/dom(2))+1)/2;
    %v0 = phi0*(exp(x-1e-6) - exp(-1e-6));
    sol0 = [c1_0, c2_0, v0];

    opts = pdeset('Eps', tol, 'PDEflag', pdeflag);

    %% Call pde15s to solve the problem.
    [t, c1, c2, v] = pde15s(pdefun, t, sol0, bc, opts);

    %% Plot the solution components.
    if plotGraphs==true
        % Plot concentration c1 (+ve charge)
        figure;
        %subplot(2,1,1);
        if intTimes==true
            k = plot(c1(:,2:end-1),'b--','LineWidth',2,'MarkerSize',20);
            hold on;
        end
        
        p = plot(c1(:,1),'k--','LineWidth',2,'MarkerSize',20);
        hold on;
        q = plot(c1(:,end),'r--','LineWidth',3,'MarkerSize',20);
        %set(gca, 'YScale', 'linear','FontSize',18);
        %title('Chebfun: Molar Concentration $c_{cations}$','interpreter','latex');
        %xlabel('Domain (x)');

        % Plot concentration c2 (-ve charge)
        %subplot(2,1,2);
        if intTimes==true
            l = plot(c2(:,2:end-1),'b','LineWidth',2,'MarkerSize',20);
        end
        %hold on
        m(1) = plot(c2(:,1),'k-','LineWidth',2,'MarkerSize',20);
        m(2) = plot(c2(:,end),'r','LineWidth',3,'MarkerSize',20);
        set(gca, 'YScale', 'linear','FontSize',18);
        title('Chebfun: Molar Concentrations $c_{i}$','interpreter','latex');
        xlabel('x [m]','interpreter','latex');
        ylabel(['Concentration ' '[$\frac{mol}{m}$]'],'interpreter','latex');
        hold off;
        if intTimes==true
            legend([k(1) p q l(1) m(1) m(2)],['$c_{1}(0<t<$' num2str(t_end) ')'],'$c_{1}(0)$',['c_{1}(' num2str(t_end) ')'] ...
                ,['$c_{2}(0<t<$' num2str(t_end) ')'],'$c_{2}(0)$',['$c_{2}$(' num2str(t_end) ')'],'Location','best','interpreter','latex');
        else
            legend([p q m(1) m(2)],'$c_{1}(0)$',['$c_{1}($' num2str(t_end) ')'] ...
                ,'$c_{2}(0)$',['$c_{2}($' num2str(t_end) ')'],'Location','best','interpreter','latex');
        end
        
        % Plot electric potential
        figure
        if intTimes==true
            k = plot(v(:,2:end-1),'b--','LineWidth',2,'MarkerSize',20);
            hold on;
        end
        p(1) = plot(v(:,1),'k-','LineWidth',2);
        hold on;
        p(2) = plot(v(:,end),'r','LineWidth',3);
        title('Chebfun: Electric Potential $\phi$','interpreter','latex');
        set(gca,'FontSize',18);
        xlabel('x [m]','interpreter','latex');
        ylabel(['Potential ' '[V]'],'interpreter','latex');
        if intTimes==true
            legend([k(1) p(1) p(2)],['0<t<' num2str(t_end)],'t=0',['t=' num2str(t_end)],'Location','best');
        else
            legend([ p(1) p(2)],'$\phi(0)$',['$\phi($' num2str(t_end) ')'],'Location','best','interpreter','latex');
        end

    end

end

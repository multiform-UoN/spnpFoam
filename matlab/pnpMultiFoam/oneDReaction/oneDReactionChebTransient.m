function [uu, vv] = oneDReactionChebTransient(tol,plotGraphs,tLinear)
    %% multiRegion.m -- an executable m-file for solving a partial differential equation

    %% Problem description.
    % Solving
    %   u_t =   D1*u" + V1u' - R1u,
    %   v_t =    D2*v" + V2v' - R2v,
    % for Neumann outer BC's and reactive interface

    %% Problem set-up
    % Create an interval of the space domain...
    % Reference length
    L0 = 1;
    L = 0;
    R = L0;
    dom = [L R];
    %...and specify a sampling of the time domain:
    t_end = 15;

    dt = 0.5;
    if(tLinear)
        t = 0:dt:t_end;
    else
        t= logspace(-3,log10(t_end),10);
    end

    %% coeffiecients
    D1 = 1;
    D2 = .1;
    V1 = 1;%0.01;
    V2 = 0;
    R1 = 0;
    R2 = 0;
    K_f = 2;
    K_r = 1;
    
    % Make the right-hand side of the PDE.
    pdefun = @(t,x,u,v) [diff(D1.*diff(u) - V1.*u)-R1.*u; ...
                        diff(D2.*diff(v) - V2.*v)-R2.*v];

    % Assign boundary conditions.
    bc = @(t,x,u,v) [  ...
                       feval(u,L)-1; feval(v,R); ...
                       feval(D2.*diff(v) - V2.*v,L) + feval(K_r.*u,R) - feval(K_f.*v,L); ...   % interface BC (flux continuity)
                       feval(D1.*diff(u) - V1.*u,R) - feval(K_f.*v,L) + feval(K_r.*u,R); ...  % 2nd interface BC (reactive flux)  
                    ];
    
    
    %[feval(D1.*diff(u) - V1.*u,L); feval(D2.*diff(v) - V2.*v,R); ... % outer BCs (no flux)
    %    feval(D1.*diff(u) + V1.*u,R) - feval(D2.*diff(v) + V2.*v,L); ...   % interface BC (flux continuity)
    %    feval(D1.*diff(u) + V1.*u,R) - feval(K_f.*v,L) + feval(K_r.*u,R)]; % 2nd interface BC (reactive flux)  

    % Construct a chebfun of the space variable on the domain,
    x = chebfun(@(x) x, dom);
    % and of the initial conditions.
    u0 = exp(-200*(x).^2);
    v0 =  0 + 0*x;
    
    sol0 = [u0, v0];
    uu0 = chebfun(@(x) u0(x+L0) ,[-L0,0]);
    vv0 = chebfun(@(x) v0(x) ,[0,L0]);

    %% Setup preferences for solving the problem.
    opts = pdeset('Eps', tol);

    %% Call pde15s to solve the problem.
    [t, u, v] = pde15s(pdefun, t, sol0, bc, opts);
    uu = chebfun(@(x) u(x+L0) ,[-L0,0]);
    vv = chebfun(@(x) v(x) ,[0,L0]);

    %% Plot the solution components.
    if(plotGraphs)
        figure;
        ku = plot(uu(:,2:3:end-1),'k--');
        hold on
        ru = plot(uu(:,end),'r--','Linewidth',2);
        bu = plot(uu(:,1),'b--','LineWidth',2);

        kv = plot(vv(:,2:3:end-1),'k');
        rv = plot(vv(:,end),'r','LineWidth',2);
        bv = plot(vv(:,1),'b','LineWidth',2);

        title('Chebfun: Multi-region scalar transport','interpreter','latex');
        xlabel('x [m]','interpreter','latex');
        legend([bu ku(1) ru bv kv(1) rv],'$c_{f}(0)$',['$c_{f}(0<t<$' num2str(t_end) '$)$'],['$c_{f}($' num2str(t_end) '$)$'],...
            '$c_{s}(0)$',['$c_{s}(0<t<$' num2str(t_end) '$)$'],['$c_{s}($' num2str(t_end) '$)$'],...
            'Location','northeast','interpreter','latex');
        ylim([0 3]);
        set(gca,'FontSize',24);
    end
    disp("At end Time, Fluid interface val == "); uu(end,end)
end
function [uu, vv] = oneDReactionChebSteady(tol,plots)

%     %% Physical constants
%     L0 = 1;
%     L = 0;
%     R = L0;
%     D1 = 10;%1;
%     D2 = 0.1;%.1;
%     V1 = 10;%100;%0.01;
%     V2 = 0;
%     R1 = 0;
%     R2 = 0;
%     K_f = 1;
%     K_r = 10;
%     Ref = 0.12533141373155; % Initial total concentration in OF
% 
%     %% Problem set-up.
%     % Define domain
%     dom = [L R];
%     
%     % Assign the differential equation to a chebop on that domain.
%     N = chebop(@(x,u,v) [diff(D1.*diff(u) - V1.*u)-R1.*u; ...
%                         diff(D2.*diff(v) - V2.*v)-R2.*v],dom);
%     % Set up the rhs of the differential equation so that N(c, p, u, v) = rhs.
%     rhs = [0;0];
%     % Assign bouyndary conditions + constraints
%     N.bc = @(x,u,v) [feval(D1.*diff(u) - V1.*u,L); feval(D2.*diff(v) - V2.*v,R); ... % outer BCs (no flux)
%                         %feval(D1.*diff(u) + V1.*u,R)...
%                          feval(D2.*diff(v) + V2.*v,L)- feval(K_r.*u,R) + feval(K_f.*v,L); ...   % interface BC (flux continuity)
%                         feval(D1.*diff(u) + V1.*u,R) - feval(K_f.*v,L) + feval(K_r.*u,R); ...  % 2nd interface BC (reactive flux)  
%                        sum(u) + sum(v) - Ref]; % Conserve total mass of system
                   
%% Physical constants
    L0 = 1;
    L = 0;
    R = L0;
    D1 = 1;
    D2 = .1;
    V1 = 1;
    V2 = 0;
    R1 = 0;
    R2 = 0;
    K_f = 2;
    K_r = 1;
    Ref = 0.1253; % Integral of initial Gaussian concentration in OF

    %% Problem set-up.
    % Define domain
    dom = [L R];
    
    % Assign the differential equation to a chebop on that domain.
    N = chebop(@(x,u,v) [diff(D1.*diff(u) - V1.*u)-R1.*u; ...
                        diff(D2.*diff(v) - V2.*v)-R2.*v],dom);
    % Set up the rhs of the differential equation so that N(c, p, u, v) = rhs.
    rhs = [0;0];
    % Assign bouyndary conditions + constraints
    N.bc = @(x,u,v) [  ...
                       feval(u,L)-1; feval(v,R); ...
                       feval(D2.*diff(v) - V2.*v,L) + feval(K_r.*u,R) - feval(K_f.*v,L); ...   % interface BC (flux continuity)
                       feval(D1.*diff(u) - V1.*u,R) - feval(K_f.*v,L) + feval(K_r.*u,R); ...  % 2nd interface BC (reactive flux)  
                    ];

   %% Setup preferences for solving the problem.
    % Create a CHEBOPPREF object for passing preferences.
    % (See 'help cheboppref' for more possible options.)
    options = cheboppref();

    % Print information to the command window while solving:
    options.display = 'iter';
    

    % Option for tolerance.
    options.bvpTol = tol;

    % Option for damping.
    options.damping = false;

    % Specify the discretization to use. Possible options are:
    %  'values' (default)
    %  'coeffs'
    %  A function handle (see 'help cheboppref' for details).
    options.discretization = 'coeffs';

    % Option for determining how long each Newton step is shown.
    if plots
        options.plotting = 0.1;
    end
    %% Solve!
    % Call solvebvp to solve the problem.
    % (With the default options, this is equivalent to u = N\rhs.)
    [u, v] = solvebvp(N, rhs, options);
    
    uu = chebfun(@(x) u(x+L0) ,[-L0,0]);
    vv = chebfun(@(x) v(x) ,[0,L0]);

    %% Plot the solution.
    if plots
        figure
        plot(uu, 'LineWidth', 2);
        hold on;
        plot(vv,'LineWidth',2);
        xlabel('x [m]','interpreter','latex'), legend('$c_{f}$','$c_{s}$','interpreter','latex');
        set(gca,'FontSize',24);
    end
end
%% chebbvp.m -- an executable m-file for solving a boundary-value problem
function [c,p,u,v] = infChannelCheb(plots,closeAll)
    if closeAll
        close all
    end
    %% Problem description.
    % Solving
    %   mu*u''=1,
    %   p'+c*v'=0,
    %   epsilon*v''=-lambda*F*c,
    %   c''+npc*c'*v'+c*v''=0,
    % for x in [0, 10^-3], subject to periodic boundary conditions.

    %% Physical constants
    epsilon = 40*8.8541878176e-12; % electrical permeability
    ele_charge = 1.6021766208e-19; % elementary charge
    boltz_const = 1.38064852e-23; % Boltzmann constant
    Ava_const = 6.022140857e+23; % Avogadro constant
    F = ele_charge*Ava_const; % Faraday constant
    T = 293; % Absolute Temperature
    mu = 10; % Dynamic viscocity
    npc = (ele_charge)/(boltz_const*T); % NP flux coeff

    %% adim version
    % npc = 1;
    %epsilon = 80;
    % F = 1;

    % System constants
    H = 1e-6; % Channel height
    lambda = 1; % Asymptotic parameter
    pref = 0; % Reference pressure
    cref = 1e-3;%1e-1;
    vref = 0.05;%1e-1;
    gradp = -mu/(H^2); % Applied external force

    %% Problem set-up.
    % Define the domain.
    dom = [0, H];

    % Assign the differential equation to a chebop on that domain.
    N = chebop(@(x,c,p,u,v) [...
                            diff(c,2)+npc.*diff(c.*diff(v)); ... % concentration
                            diff(p,1)+(F.*c.*diff(v)); ...       % pressure
                            mu.*diff(u,2)-gradp; ...             % velocity
                            -epsilon.*diff(v,2)-lambda.*F.*c ...  % potential
                            ], dom);

    % Set up the rhs of the differential equation so that N(c, p, u, v) = rhs.
    rhs = [0;0;0;0];

    % Assign boundary conditions to the chebop.
    %N.lbc = @(c,p,u,v) [(c); p; u; v];   % here you can try with neumann for c
    %N.rbc = @(c,p,u,v) [c-cref; p; u; v-vref];
    
    N.bc = @(x,c,p,u,v) [ ... %feval(c,0); feval(c,H)-cref;  ...    % concentration (dirichlet)
        feval(diff(c) + npc.*c.*diff(v) , 0); feval(diff(c) + npc.*c.*diff(v), H); mean(c)- cref; ... %concentration (noflux)
        feval(p,H)-pref; ...                      % pressure
        feval(u,0); feval(u,H); ...               % velocity
        feval(v,0); feval(v,H)-vref ];            % potential
    
    % --------- Dirichlet Concentration BCs ----------
    %N.bc = @(x,c,p,u,v) [feval(c,0); feval(p,H)-pref; feval(u,0); feval(v,0); ...
    %    feval(c,H)-cref;  feval(u,H); feval(v,H)-vref];
    
    % --------- No-flux Concentration BCs ------------
   % N.bc = @(x,c,p,u,v) [feval(diff(c) + npc.*c.*diff(v),0); feval(p,H)-pref; feval(u,0); feval(v,0); ...
   %    feval(diff(c) + npc.*c.*diff(v),H); feval(u,H); feval(v,H)-vref];

    %% Setup preferences for solving the problem.
    % Create a CHEBOPPREF object for passing preferences.
    % (See 'help cheboppref' for more possible options.)
    options = cheboppref();

    % Print information to the command window while solving:
    options.display = 'iter';
    

    % Option for tolerance.
    options.bvpTol = 5e-13;

    % Option for damping.
    options.damping = false;

    % Specify the discretization to use. Possible options are:
    %  'values' (default)
    %  'coeffs'
    %  A function handle (see 'help cheboppref' for details).
    options.discretization = 'values';

    % Option for determining how long each Newton step is shown.
    if plots
        options.plotting = 0.1;
    end
    %% Solve!
    % Call solvebvp to solve the problem.
    % (With the default options, this is equivalent to u = N\rhs.)
    [c, p, u, v] = solvebvp(N, rhs, options);

    %% Plot the solution.
    if plots
        figure
        plot([c, p, u, v], 'LineWidth', 2)
        title(['Final solution lambda=',num2str(lambda)]), xlabel('x'), legend('c','p','u','v')
        figure
        plot(c/norm(c))
        hold on
        plot(p/norm(p))
        plot(u/norm(u))
        plot(v/norm(v))
        legend(['c/cnorm';'p/pnorm';'u/unorm';'v/vnorm'])
    end
end

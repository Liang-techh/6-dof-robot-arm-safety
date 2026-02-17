function visualize_6dof_arm()
    % VISUALIZE_6DOF_ARM Interactive 3D visualization with Dynamic Constraints
    % Features: Manual Control, Infinite Simulation, Perturbation Injection
    
    clear; close all; clc;

    % Load results for constraints
    data = [];
    try
        data = load('result_6dof_arm.mat');
    catch
        warning('result_6dof_arm.mat not found. Using default parameters.');
    end

    % Safe set parameters
    eta = 0.8; 
    if isfield(data, 'eta_star'), eta = data.eta_star; end
    
    % Figure Setup
    hFig = figure('Name', '6-DOF Arm Visualization System', 'NumberTitle', 'off', ...
                  'Position', [50, 50, 1400, 800], 'Color', 'w');
    
    % State Variables
    sim_running = false;
    x_state = zeros(12,1); % [q; dq]
    
    % --- Layout ---
    % Left: 3D Arm
    hAx3D = axes('Parent', hFig, 'Position', [0.05 0.1 0.45 0.85]);
    grid(hAx3D, 'on'); hold(hAx3D, 'on'); axis(hAx3D, 'equal');
    xlabel(hAx3D, 'X'); ylabel(hAx3D, 'Y'); zlabel(hAx3D, 'Z');
    view(hAx3D, 45, 30);
    xlim(hAx3D, [-0.6 0.6]); ylim(hAx3D, [-0.6 0.6]); zlim(hAx3D, [0 1.0]);
    title(hAx3D, '3D Configuration');

    % Right: Constraints
    hAxC1 = axes('Parent', hFig, 'Position', [0.6 0.7 0.35 0.22]);
    hAxC2 = axes('Parent', hFig, 'Position', [0.6 0.4 0.35 0.22]);
    hAxC3 = axes('Parent', hFig, 'Position', [0.6 0.1 0.35 0.22]);
    setupConstraintAxis(hAxC1, 'q1', 'q2');
    setupConstraintAxis(hAxC2, 'q3', 'q4');
    setupConstraintAxis(hAxC3, 'q5', 'q6');

    % Draw Static Boundaries
    drawSafetySet(hAxC1, eta, 1.5);
    drawSafetySet(hAxC2, eta, 1.5);
    drawSafetySet(hAxC3, eta, 1.5);

    % --- Plot Objects ---
    hLine = plot3(hAx3D, 0, 0, 0, 'b-o', 'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    hEffector = plot3(hAx3D, 0, 0, 0, 'rx', 'LineWidth', 2, 'MarkerSize', 10);
    hSafeDots = [
        plot(hAxC1, 0, 0, 'ro', 'MarkerFaceColor', 'r');
        plot(hAxC2, 0, 0, 'ro', 'MarkerFaceColor', 'r');
        plot(hAxC3, 0, 0, 'ro', 'MarkerFaceColor', 'r')
    ];
    hTrails = [
        plot(hAxC1, 0, 0, 'r-', 'LineWidth', 1);
        plot(hAxC2, 0, 0, 'r-', 'LineWidth', 1);
        plot(hAxC3, 0, 0, 'r-', 'LineWidth', 1)
    ];

    % --- Controls ---
    hBtnSim = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
        'String', 'Start Simulation (Infinite)', ...
        'Position', [50, 50, 150, 40], ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'Callback', @toggleSimulation);

    hBtnPerturb = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
        'String', 'Inject Disturbance', ...
        'Position', [220, 50, 150, 40], ...
        'FontSize', 10, ...
        'Callback', @injectDisturbance, ...
        'Enable', 'off');

    % Sliders
    hSld = gobjects(6,1); hLbl = gobjects(6,1);
    for i = 1:6
        hLbl(i) = uicontrol('Parent', hFig, 'Style', 'text', ...
            'Position', [50, (i-1)*35 + 110, 60, 20], ...
            'String', sprintf('q%d: 0.00', i), 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
        hSld(i) = uicontrol('Parent', hFig, 'Style', 'slider', ...
            'Position', [120, (i-1)*35 + 110, 300, 20], ...
            'Min', -pi, 'Max', pi, 'Value', 0, ...
            'Callback', @sliderCallback);
    end

    % DH Params & Dynamics
    L_BASE = 0.10; D_BASE = 0.08; L_ARM = 0.21; L_FOREARM = 0.19; D_ELBOW = 0.05; L_WRIST = 0.07;
    DH_params = [
         0,          L_BASE,    D_BASE,   -pi/2;
        -pi/2,       0,         L_ARM,     0;
         pi/2,       D_ELBOW,   0,         pi/2;
         0,          L_FOREARM, 0,        -pi/2;
         0,          0,         0,         pi/2;
         0,          L_WRIST,   0,         0
    ];

    a_coef = [1.80; 2.17; 2.86; 3.75; 5.80; 8.60];
    c_coef = [1.30; 1.83; 2.71; 4.00; 6.50; 10.00];
    gw = [1.0; 0.5; 0.3; 0.2; 0.1; 0.05];
    coupling = 0.05;

    if isfield(data, 'a_coef'), a_coef = data.a_coef; end
    if isfield(data, 'c_coef'), c_coef = data.c_coef; end
    if isfield(data, 'gw'), gw = data.gw; end
    if isfield(data, 'coupling'), coupling = data.coupling; end

    % Dynamics Matrix
    A_mat = diag(-a_coef);
    A_mat(1,2)=coupling; A_mat(2,1)=coupling; A_mat(2,3)=coupling; A_mat(3,2)=coupling; 
    A_mat(3,4)=coupling; A_mat(4,3)=coupling; A_mat(4,5)=coupling; A_mat(5,4)=coupling; 
    A_mat(5,6)=coupling; A_mat(6,5)=coupling;
    C_mat = diag(-c_coef);
    odefun = @(tv, xv) [xv(7:12); A_mat*xv(1:6) + C_mat*xv(7:12)]; % + disturbance handled manually

    % Initialize
    updateVisuals(zeros(6,1));

    % --- Callbacks ---

    function sliderCallback(~, ~)
        if sim_running, return; end % Ignore if simulating
        q_val = zeros(6,1);
        for k=1:6
            q_val(k) = hSld(k).Value;
            hLbl(k).String = sprintf('q%d: %.2f', k, q_val(k));
        end
        x_state(1:6) = q_val;
        x_state(7:12) = 0; % Stop moving if slider touched
        updateVisuals(q_val);
    end

    function toggleSimulation(~, ~)
        sim_running = ~sim_running;
        if sim_running
            hBtnSim.String = 'STOP Simulation';
            hBtnSim.BackgroundColor = [1 0.7 0.7];
            hBtnPerturb.Enable = 'on';
            disableSliders();
            
            % If currently at zero, start random
            if norm(x_state) < 1e-3
                randomStart(); 
            end
            
            runSimulationLoop();
        else
            hBtnSim.String = 'Start Simulation (Infinite)';
            hBtnSim.BackgroundColor = [0.94 0.94 0.94];
            hBtnPerturb.Enable = 'off';
            enableSliders();
        end
    end

    function injectDisturbance(~, ~)
        % Kick velocity states
        kick = (rand(6,1)-0.5)*2.0; % Random kick
        x_state(7:12) = x_state(7:12) + kick;
    end

    function randomStart()
        R_init = 0.15; if isfield(data, 'R_init'), R_init = data.R_init; end
        dir = randn(6,1); dir = dir/norm(dir);
        x_state = [dir * R_init * 0.9; zeros(6,1)];
    end

    function runSimulationLoop()
        dt = 0.05;
        trail_buf = zeros(12, 100); % Keep last 100 points
        trail_idx = 1;
        
        while sim_running && ishandle(hFig)
            % Integrate one step (RK4 or ODE45 step)
            % Simple Euler for responsiveness: x = x + dx*dt
            % Or better: ode45 for one small chunk
            [~, Y] = ode45(odefun, [0 dt], x_state);
            x_state = Y(end, :)';
            
            % Update Visuals
            q_curr = x_state(1:6);
            updateVisuals(q_curr);
            
            % Update Trails
            trail_buf(:, trail_idx) = x_state;
            trail_idx = mod(trail_idx, 100) + 1;
            
            % Extract valid trail points for plotting
            % We just plot the whole buffer but it might have jumps if not careful.
            % Better: just plot what we have so far (up to 100 pts).
            % A moving buffer is nicer.
            
            % Quick hack for visual trail: just last N points
            
            set(hTrails(1), 'XData', trail_buf(1,:), 'YData', trail_buf(2,:));
            set(hTrails(2), 'XData', trail_buf(3,:), 'YData', trail_buf(4,:));
            set(hTrails(3), 'XData', trail_buf(5,:), 'YData', trail_buf(6,:));
             
            % Update Sliders to reflect simulation
            for k=1:6
                hSld(k).Value = max(min(q_curr(k), pi), -pi);
                hLbl(k).String = sprintf('q%d: %.2f', k, q_curr(k));
            end

            drawnow;
        end
        
        if ishandle(hFig)
            hBtnSim.String = 'Start Simulation (Infinite)';
            hBtnSim.BackgroundColor = [0.94 0.94 0.94];
            hBtnPerturb.Enable = 'off';
            enableSliders();
        end
    end

    % --- Helpers ---
    function enableSliders()
        for k=1:6, hSld(k).Enable = 'on'; end
    end

    function disableSliders()
        for k=1:6, hSld(k).Enable = 'off'; end
    end

    function updateVisuals(q)
        [pts, ~] = forward_kinematics_6dof(q, DH_params);
        set(hLine, 'XData', pts(1,:), 'YData', pts(2,:), 'ZData', pts(3,:));
        set(hEffector, 'XData', pts(1,end), 'YData', pts(2,end), 'ZData', pts(3,end));
    end

    function setupConstraintAxis(hAx, xlab, ylab)
        grid(hAx, 'on'); hold(hAx, 'on'); axis(hAx, 'equal');
        xlabel(hAx, [xlab ' (rad)']); ylabel(hAx, [ylab ' (rad)']);
        xlim(hAx, [-1 1]); ylim(hAx, [-1 1]);
    end

    function drawSafetySet(hAx, eta, p_coeff)
        r = sqrt(eta/p_coeff);
        th = linspace(0, 2*pi, 100);
        xc = r*cos(th); yc = r*sin(th);
        plot(hAx, xc, yc, 'k-', 'LineWidth', 2);
        fill(hAx, xc, yc, [0.9 1 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end

end

function [points, T_end] = forward_kinematics_6dof(q, dh)
    num_joints = 6;
    T = eye(4);
    points = zeros(3, num_joints + 1);
    points(:,1) = [0; 0; 0];
    for i = 1:num_joints
        theta = q(i) + dh(i,1);
        d = dh(i,2); a = dh(i,3); alpha = dh(i,4);
        ct = cos(theta); st = sin(theta); ca = cos(alpha); sa = sin(alpha);
        T_i = [ct -st*ca st*sa a*ct; st ct*ca -ct*sa a*st; 0 sa ca d; 0 0 0 1];
        T = T * T_i;
        points(:, i+1) = T(1:3, 4);
    end
    T_end = T;
end

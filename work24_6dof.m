clearvars; close all; clc;
yalmip('clear');

%% ==== 6-DOF Robotic Arm SOS Reachability ====
% 12 states: x = [q1..q6, dq1..dq6]
% DH parameters from 6dof_kinematic.cpp
% Linearized dynamics for SOS tractability (degV=2)

disp('==== 6-DOF Arm: Parameters ====');

% DH-based link parameters (normalized to ~O(1))
% From 6dof_kinematic.cpp: L_BASE, D_BASE, L_ARM, L_FOREARM, D_ELBOW, L_WRIST
L = [0.10; 0.08; 0.21; 0.19; 0.05; 0.07];  % link lengths (m)

% Masses (decreasing from base to tip)
m = [1.0; 0.8; 0.6; 0.4; 0.3; 0.15];

% Inertia (I ~ m*L^2, normalized)
I_val = [1.0; 0.6; 0.35; 0.2; 0.1; 0.05];

% Gravity torque (mgl at CoM, decreasing distally)
mgl = [0.8; 0.5; 0.3; 0.15; 0.08; 0.03];

% Friction
b_fr = [0.5; 0.4; 0.35; 0.3; 0.25; 0.2];

% PD gains (strong for stability)
Kp = [1.0; 0.8; 0.7; 0.6; 0.5; 0.4];
Kd = [0.8; 0.7; 0.6; 0.5; 0.4; 0.3];

% Derived: linear coefficients per joint
a_coef = (Kp + mgl) ./ I_val;   % restoring
c_coef = (b_fr + Kd) ./ I_val;  % damping
gw = [1.0; 0.5; 0.3; 0.2; 0.1; 0.05] ./ I_val;  % disturbance gain

% Coupling (nearest-neighbor, small)
coupling = 0.05;

for j=1:6
    fprintf('  Joint %d: a=%.2f, c=%.2f, gw=%.3f\n', j, a_coef(j), c_coef(j), gw(j));
end

%% ==== Dynamics ====
disp('==== Defining 12-state Linear Dynamics ====');

T = 1; t0 = 0;
R_init = 0.15;
degV = 2;
degS = 0;   % Constant multipliers for 12-state tractability

sdpvar t q1 q2 q3 q4 q5 q6 dq1 dq2 dq3 dq4 dq5 dq6 w real
x = [q1;q2;q3;q4;q5;q6;dq1;dq2;dq3;dq4;dq5;dq6];
q = [q1;q2;q3;q4;q5;q6];
dq = [dq1;dq2;dq3;dq4;dq5;dq6];

% f = [dq; ddq] — linearized E-L with nearest-neighbor coupling
% Define each ddq explicitly (cannot use zeros() with sdpvar)
ddq1 = -a_coef(1)*q1 - c_coef(1)*dq1 + coupling*q2 + gw(1)*w;
ddq2 = -a_coef(2)*q2 - c_coef(2)*dq2 + coupling*q1 + coupling*q3 + gw(2)*w;
ddq3 = -a_coef(3)*q3 - c_coef(3)*dq3 + coupling*q2 + coupling*q4 + gw(3)*w;
ddq4 = -a_coef(4)*q4 - c_coef(4)*dq4 + coupling*q3 + coupling*q5 + gw(4)*w;
ddq5 = -a_coef(5)*q5 - c_coef(5)*dq5 + coupling*q4 + coupling*q6 + gw(5)*w;
ddq6 = -a_coef(6)*q6 - c_coef(6)*dq6 + coupling*q5 + gw(6)*w;

f = [dq1;dq2;dq3;dq4;dq5;dq6; ddq1;ddq2;ddq3;ddq4;ddq5;ddq6];
disp('  12-state linearized dynamics defined.');

%% ==== SOS Setup ====
disp('==== SOS Variable Setup ====');

[V, cV] = polynomial([t; x], degV);

g_time = (t - t0)*(T - t);

% Level-set polynomials
p_poly = 1.5*(q'*q) + 0.8*(dq'*dq);
q_poly = 3.0*(q'*q) + 2.0*(dq'*dq);

r0poly = x'*x - R_init^2;

% SOS multipliers (constant = scalar, degS=0)
[s1,cS1] = polynomial([x; w; t], degS);
[s2,cS2] = polynomial([x; w; t], degS);
[s4,cS4] = polynomial(x, degS);
[s5,cS5] = polynomial(x, degS);
[s6,cS6] = polynomial([x; t], degS);
[s7,cS7] = polynomial([x; t], degS);

dVdt = jacobian(V, t);
gradV = jacobian(V, x);
Vdot = dVdt + gradV * f;
Vt0 = replace(V, t, t0);

R = 1;
h_t = t^2 / T^2;

allCoeff = [cV; cS1; cS2; cS4; cS5; cS6; cS7];
options = sdpsettings('solver', 'mosek', 'verbose', 1, ...
                      'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 600);

fprintf('  V coefficients: %d\n', length(cV));
fprintf('  Total decision vars: %d\n', length(allCoeff));

%% ==== Bisection on eta ====
eta_min = 0.01; eta_max = 3.0; eta_tol = 0.05;
max_iter = 15; eta_star = NaN;
disp('==== Bisection on eta ====');

for iter = 1:max_iter
    eta_try = (eta_min + eta_max) / 2;
    cons = [];
    
    diss = -(Vdot - w^2) + (p_poly - eta_try)*s1 - s2*g_time;
    cons = [cons, sos(diss)];
    cons = [cons, sos(-Vt0 + s4*r0poly)];
    cons = [cons, sos(-(p_poly - eta_try)*s6 + V - R^2*h_t - s7*g_time)];
    cons = [cons, sos(s1), sos(s2), sos(s4), sos(s6), sos(s7)];
    
    try
        sol = solvesos(cons, [], options, allCoeff);
        is_feasible = (sol.problem == 0);
    catch ME
        disp(['  Error: ', ME.message]);
        is_feasible = false;
    end
    
    if is_feasible
        fprintf('  Iter %d: eta=%.4f -> FEASIBLE\n', iter, eta_try);
        eta_max = eta_try; eta_star = eta_try;
    else
        fprintf('  Iter %d: eta=%.4f -> infeasible\n', iter, eta_try);
        eta_min = eta_try;
    end
    if abs(eta_max - eta_min) < eta_tol, break; end
end

% Fallback
if isnan(eta_star)
    warning('Trying R_init=0.08...');
    R_init = 0.08;
    r0poly = x'*x - R_init^2;
    Vt0 = replace(V, t, t0);
    eta_min = 0.01; eta_max = 2.0;
    for iter = 1:max_iter
        eta_try = (eta_min + eta_max) / 2;
        cons = [];
        diss = -(Vdot - w^2) + (p_poly - eta_try)*s1 - s2*g_time;
        cons = [cons, sos(diss)];
        cons = [cons, sos(-Vt0 + s4*r0poly)];
        cons = [cons, sos(-(p_poly - eta_try)*s6 + V - R^2*h_t - s7*g_time)];
        cons = [cons, sos(s1), sos(s2), sos(s4), sos(s6), sos(s7)];
        try
            sol = solvesos(cons, [], options, allCoeff);
            is_feasible = (sol.problem == 0);
        catch, is_feasible = false; end
        if is_feasible
            fprintf('  Fallback %d: eta=%.4f -> FEASIBLE\n', iter, eta_try);
            eta_max = eta_try; eta_star = eta_try;
        else
            fprintf('  Fallback %d: eta=%.4f -> infeasible\n', iter, eta_try);
            eta_min = eta_try;
        end
        if abs(eta_max - eta_min) < eta_tol, break; end
    end
end

if isnan(eta_star)
    error('6-DOF: Could not find feasible eta.');
end
disp(['  --> 6-DOF eta* = ', num2str(eta_star)]);

%% ==== Bisection on alpha ====
alpha_min = 0.01; alpha_max = 8.0; alpha_tol = 0.05;
alpha_star = NaN;
disp('==== Bisection on alpha ====');

for iter = 1:max_iter
    alpha_try = (alpha_min + alpha_max) / 2;
    cons = [];
    
    diss = -(Vdot - w^2) + (p_poly - eta_star)*s1 - s2*g_time;
    cons = [cons, sos(diss)];
    cons = [cons, sos(-Vt0 + s4*r0poly)];
    cons = [cons, sos(-(p_poly - eta_star)*s6 + V - R^2*h_t - s7*g_time)];
    
    VT = replace(V, t, T);
    cons = [cons, sos(-(q_poly - alpha_try)*s5 + VT - R^2)];
    
    cons = [cons, sos(s1), sos(s2), sos(s4), sos(s5), sos(s6), sos(s7)];
    
    try
        sol = solvesos(cons, [], options, allCoeff);
        is_feasible = (sol.problem == 0);
    catch, is_feasible = false; end
    
    if is_feasible
        fprintf('  Iter %d: alpha=%.4f -> FEASIBLE\n', iter, alpha_try);
        alpha_max = alpha_try; alpha_star = alpha_try;
    else
        fprintf('  Iter %d: alpha=%.4f -> infeasible\n', iter, alpha_try);
        alpha_min = alpha_try;
    end
    if abs(alpha_max - alpha_min) < alpha_tol, break; end
end

if isnan(alpha_star)
    warning('Alpha failed, using 2*eta.');
    alpha_star = eta_star * 2;
end
disp(['  --> 6-DOF alpha* = ', num2str(alpha_star)]);

%% ==== Extract V ====
disp('==== Extracting V(t,x) ====');

syms ts qs1 qs2 qs3 qs4 qs5 qs6 dqs1 dqs2 dqs3 dqs4 dqs5 dqs6 real
x_sym = [qs1;qs2;qs3;qs4;qs5;qs6;dqs1;dqs2;dqs3;dqs4;dqs5;dqs6];
monoms_sym = monolist([ts; x_sym], degV);
cV_val = value(cV);
Vsym = cV_val' * monoms_sym;
disp('  V extracted.');

clear cV cS1 cS2 cS4 cS5 cS6 cS7 s1 s2 s4 s5 s6 s7
clear sol cons allCoeff V Vt0 VT Vdot

save('result_6dof_arm.mat', 'alpha_star', 'eta_star', 'Vsym', 'R_init', 'T', 'degV');

%% ==== ODE45 Simulation ====
disp('==== ODE45 Simulations ====');

A_mat = diag(-a_coef);
A_mat(1,2) = coupling; A_mat(2,1) = coupling;
A_mat(2,3) = coupling; A_mat(3,2) = coupling;
A_mat(3,4) = coupling; A_mat(4,3) = coupling;
A_mat(4,5) = coupling; A_mat(5,4) = coupling;
A_mat(5,6) = coupling; A_mat(6,5) = coupling;

C_mat = diag(-c_coef);

odefun = @(tv, xv, wv) [
    xv(7:12);
    A_mat*xv(1:6) + C_mat*xv(7:12) + gw*wv
];

num_sim = 500;
x_T = zeros(num_sim, 12);
rng(42);
c_bound = sqrt(3) * 0.2;

for k = 1:num_sim
    dir = randn(12,1); dir = dir/norm(dir);
    r = R_init * rand^(1/12) * 0.95;
    x0 = r * dir;
    c_k = (rand*2-1)*c_bound;
    [~, sol_k] = ode45(@(tv,xv) odefun(tv,xv,c_k*tv), [t0 T], x0);
    x_T(k,:) = sol_k(end,:);
end
disp(['  ', num2str(num_sim), ' simulations done.']);

%% ==== Visualization ====
disp('==== Visualization ====');

figure('Position', [50 50 1600 500], 'Color', 'w');

% --- q1 vs q2 slice ---
subplot(1,3,1); hold on; grid on;
xr = linspace(-0.8, 0.8, 200);
[X1,X2] = meshgrid(xr, xr);

VT_12 = subs(Vsym, ts, T);
for j = 3:6, VT_12 = subs(VT_12, x_sym(j), 0); end
for j = 7:12, VT_12 = subs(VT_12, x_sym(j), 0); end
VT12_fun = matlabFunction(VT_12, 'Vars', {qs1, qs2});
VT12_v = arrayfun(VT12_fun, X1, X2);

p12 = @(a,b) 1.5*(a.^2+b.^2);
contour(X1,X2,p12(X1,X2),[eta_star eta_star],'Color',[0.93 0.69 0.13],'LineWidth',2);
contour(X1,X2,3*(X1.^2+X2.^2),[alpha_star alpha_star],'b--','LineWidth',2);
contour(X1,X2,VT12_v,[R^2 R^2],'k-','LineWidth',2.5);
th=linspace(0,2*pi,200);
plot(R_init*cos(th),R_init*sin(th),'r:','LineWidth',1.5);
plot(x_T(:,1),x_T(:,2),'m.','MarkerSize',3);

xlabel('$q_1$','Interpreter','latex'); ylabel('$q_2$','Interpreter','latex');
title('$q_1$-$q_2$','Interpreter','latex','FontSize',13);
axis equal; xlim([-0.8 0.8]); ylim([-0.8 0.8]);

% --- q1 vs q4 slice ---
subplot(1,3,2); hold on; grid on;
VT_14 = subs(Vsym, ts, T);
for j = [2,3,5,6], VT_14 = subs(VT_14, x_sym(j), 0); end
for j = 7:12, VT_14 = subs(VT_14, x_sym(j), 0); end
VT14_fun = matlabFunction(VT_14, 'Vars', {qs1, qs4});
VT14_v = arrayfun(VT14_fun, X1, X2);

contour(X1,X2,p12(X1,X2),[eta_star eta_star],'Color',[0.93 0.69 0.13],'LineWidth',2);
contour(X1,X2,3*(X1.^2+X2.^2),[alpha_star alpha_star],'b--','LineWidth',2);
contour(X1,X2,VT14_v,[R^2 R^2],'k-','LineWidth',2.5);
plot(R_init*cos(th),R_init*sin(th),'r:','LineWidth',1.5);
plot(x_T(:,1),x_T(:,4),'m.','MarkerSize',3);

xlabel('$q_1$','Interpreter','latex'); ylabel('$q_4$','Interpreter','latex');
title('$q_1$-$q_4$','Interpreter','latex','FontSize',13);
axis equal; xlim([-0.8 0.8]); ylim([-0.8 0.8]);

% --- q1 vs dq1 slice ---
subplot(1,3,3); hold on; grid on;
VT_1d1 = subs(Vsym, ts, T);
for j = 2:6, VT_1d1 = subs(VT_1d1, x_sym(j), 0); end
for j = 8:12, VT_1d1 = subs(VT_1d1, x_sym(j), 0); end
VT1d1_fun = matlabFunction(VT_1d1, 'Vars', {qs1, dqs1});
VT1d1_v = arrayfun(VT1d1_fun, X1, X2);

p1d1 = @(a,b) 1.5*a.^2 + 0.8*b.^2;
contour(X1,X2,p1d1(X1,X2),[eta_star eta_star],'Color',[0.93 0.69 0.13],'LineWidth',2);
contour(X1,X2,3*X1.^2+2*X2.^2,[alpha_star alpha_star],'b--','LineWidth',2);
contour(X1,X2,VT1d1_v,[R^2 R^2],'k-','LineWidth',2.5);
plot(R_init*cos(th),R_init*sin(th),'r:','LineWidth',1.5);
plot(x_T(:,1),x_T(:,7),'m.','MarkerSize',3);

xlabel('$q_1$','Interpreter','latex'); ylabel('$\dot{q}_1$','Interpreter','latex');
title('$q_1$-$\dot{q}_1$','Interpreter','latex','FontSize',13);
axis equal; xlim([-0.8 0.8]); ylim([-0.8 0.8]);

savefig('result_6dof_arm.fig');
print('result_6dof_arm', '-dpng', '-r150');
disp('Figures saved.');

%% ==== Summary ====
fprintf('\n============================================\n');
fprintf('  6-DOF Robotic Arm SOS Reachability\n');
fprintf('============================================\n');
fprintf('  eta*   = %.4f\n', eta_star);
fprintf('  alpha* = %.4f\n', alpha_star);
fprintf('  R_init = %.3f\n', R_init);
fprintf('  degV=%d, degS=%d, states=12\n', degV, degS);
fprintf('\nProgression:\n');
fprintf('  1-DOF:  2 states, degV=8\n');
fprintf('  2-DOF:  4 states, degV=4\n');
fprintf('  3-DOF:  6 states, degV=2\n');
fprintf('  6-DOF: 12 states, degV=2 (linearized)\n');
fprintf('============================================\n');
disp('==== ALL DONE ====');

clc;clear all; close all
import casadi.*
x= SX.sym('w'); %decision Variables
obj = exp(0.2*x)*sin(x); %Calculate Obj
g = []; % Optimization Constraints - Empty(Unconstrained)
p = []; %Optimization Problem Parameters - Empty(No parameters used here)

opt_variables = x; % Single Decision Variable
nlp_prob = struct('f', obj, 'x', opt_variables, 'g', g, 'p', p);
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;%0,3
opts.print_time = 0;%0,1
opts.ipopt.acceptable_tol = 1e-8; %Optimality Convergence Tolerance 
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;
args.lbx = 0;
args.ubx = 4*pi;
args.lbg = -inf;
args.ubg = inf;

args.p = []; %There are no parameters in this optimization problem
args.x0 = 4; % Initialization of the Optimization 

sol = solver('x0',args.x0,'lbx',args.lbx, 'ubx',args.ubx, 'lbg',args.lbg, 'ubg',args.ubg, 'p',args.p);

sol_x = full(sol.x);
min_value = full(sol.f);
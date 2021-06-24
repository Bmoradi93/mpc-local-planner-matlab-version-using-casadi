clc; clear all; close all
import casadi.*

% mypi = raspi();
% srv = servo (mypi, 26);
% writePosition(srv, 90);
% 
% configurePin(mypi, 12, 'PWM');
% configurePin(mypi, 21, 'PWM');
% 
% writePWMFrequency(mypi, 12, 2000);
% writePWMFrequency(mypi, 21, 2000);
% 
% writePWMVoltage(mypi, 12, 0);
% writePWMVoltage(mypi, 21, 0);

T = 0.2; % sampling time [s]
N = 20; % prediction horizon
rob_diam = 0.65;

v_max = 0.9; 
v_min = -v_max;

omega_max = pi/3; 
omega_min = -omega_max;

x = SX.sym('x'); 
y = SX.sym('y'); 
theta = SX.sym('theta');
v = SX.sym('v'); 
omega = SX.sym('omega');

states = [x;y;theta]; 
n_states = length(states);

controls = [v;omega]; 
n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs});
U = SX.sym('U',n_controls,N);
P = SX.sym('P',n_states + n_states);
X = SX.sym('X',n_states,(N+1));

% compute solution symbolically
X(:,1) = P(1:3); % initial state
for k = 1:N
    st = X(:,k);  
    con = U(:,k);
    f_value  = f(st,con);
    st_next  = st+ (T*f_value);
    X(:,k+1) = st_next;
end
% this function to get the optimal trajectory knowing the optimal solution
ff=Function('ff',{U,P},{X});

obj = 0; % Objective function
Q = zeros(3,3); 
Q(1,1) = 1;
Q(2,2) = 5;
Q(3,3) = 0.1; % weighing matrices (states)

R = zeros(2,2); 
R(1,1) = 0.5; 
R(2,2) = 0.05; % weighing matrices (controls)
% compute objective
for k=1:N
    st = X(:,k);  
    con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
end

% compute constraints
g = [];  % constraints vector
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; X(1,k)];   %state x
    g = [g ; X(2,k)];   %state y
end

% make the decision variables one column vector
OPT_variables = reshape(U,2*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;
% inequality constraints (state constraints)
args.lbg = -2;  % lower bound of the states x and y
args.ubg = 2;   % upper bound of the states x and y 

% input constraints
args.lbx(1:2:2*N-1,1) = v_min; 
args.lbx(2:2:2*N,1)   = omega_min;
args.ubx(1:2:2*N-1,1) = v_max; 
args.ubx(2:2:2*N,1)   = omega_max;


%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0 ; 0 ; 0.0];    % initial condition.
xs = [1.5; 1.5 ; 0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);  % two control inputs 

sim_tim = 20; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];


% the main simulaton loop... it works as long as the error is greater
% than 10^-2 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 1e-2 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    args.x0 = reshape(u0',2*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',2,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:3,mpciter+1)= full(ff_value)';
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Deploy the Signal--------------------------
%     writePosition(srv, 90-(u(1,2)*180/pi));
%     writePWMVoltage(mypi, 12, abs(1.66*u(1,1)));
%     writePWMVoltage(mypi, 21, abs(1.66*u(1,1)));
    %--------------------------------------------
    [t0, x0, u0] = shift(T, t0, x0, u,f); % get the initialization of the next optimization step
    
    xx(:,mpciter+2) = x0;
    mpciter;
    mpciter = mpciter + 1;
end

% Results
figure(1);
plot(xx(1,:), xx(2,:), 'LineWidth', 2, 'color', 'red');
grid on
box on
title('The Simulated Trajectory for the Robot');
xlabel('X(m)');
ylabel('Y(m)');

figure(2);
subplot(211);
plot(u_cl(:,1), 'LineWidth', 2, 'color', 'blue');
grid on
box on
title('The Input Speed of the Robot');
xlabel('T(s)');
ylabel('V(m/s)');

% figure(3);
subplot(212);
plot(xx(3,:), 'LineWidth', 2, 'color', 'black');
grid on
box on
title('The Input Angle of the Robot');
xlabel('T(s)');
ylabel('Theta(m/s)');
% main_loop_time = toc(main_loop);
% ss_error = norm((x0-xs),2);
% average_mpc_time = main_loop_time/(mpciter+1);
% 
% Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam) % a drawing function

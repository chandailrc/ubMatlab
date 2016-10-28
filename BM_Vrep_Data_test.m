% CONSTANTS
% datafile = csvread('dataForMatlab_Moving.csv');
datafile = csvread('dataForMatlab_static.csv');
dt = 0.15;
m = 941.20; 	% Mass of vehicle in kgg
h = 0.84; 	% m
tf = 1.28;	% m
tr = tf;
a = 0.968;	% m
b = 1.392;	% m
g = 9.81;
NN = 1000;
Iz = 1125.67;
Ca_FL = 40000.0; %200520.0    # cornering stiffness / lateral stiffness
Ca_FR = 40000.0;
Ca_RL = 40000.0; %76426.0
Ca_RR = 40000.0;
Ck_FL = 150000.0;  %30.0        # braking stiffness / longitudinal stiffness
Ck_FR = 150000.0;  %30.0
Ck_RL = 150000.0;  %24.0
Ck_RR = 150000.0;  %24.0
wheelr = (0.4/2.0) + 0.8;
% Number of States
n_states=7;
n_meas_states = 7;

% Init state

s_state=[...
        datafile(32,2);...    %         datafile(32,2);...         % 1. x
        datafile(32,3);...    %         datafile(32,3);...         % 2. y
        datafile(32,11);...   %         datafile(32,11);...        % 3. vxz
        datafile(32,12);...   %         datafile(32,12);...        % 4. vy
        datafile(32,14);...   %         datafile(32,14);...        % 5. ax
        datafile(32,15);...   %         datafile(32,15);...        % 6. ay
        datafile(32,5);...    %         datafile(32,5);...         % 7. yaw
        ];

    
% % Q AND R
% q_proc= 10;                                 %std of process 
% r_meas= 0;                                 %std of measurement
% Q_proc=(q_proc^2)*eye(n_states);              % covariance of process
% % Q_proc = diag([]);
% R_meas=r_meas^2*eye(n_meas_states);

Q_proc = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.001]);
% R_meas=r_meas^2*eye(n_meas_states);
R_meas = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.001]);

    
x=s_state+q_proc*0.01*randn(n_states,1);
P = 10*eye(n_states);
x_store = zeros(NN,n_states);
x1_store = zeros(NN,n_states);%estmate
meas_store = zeros(n_meas_states,NN);             %act





for k=32: 1031
    meas = [datafile(k, 2); datafile(k,3); datafile(k,11); datafile(k,12); datafile(k,14); datafile(k,15); datafile(k,5)];
    delta = datafile(k,17);
    awsFL = datafile(k,18);
    awsFR = datafile(k,19);
    awsRL = datafile(k,20);
    awsRR = datafile(k,21);
    
    [x, P, x1] = BM_ekf(x, P, Q_proc, R_meas, meas, delta, awsFL, awsFR, awsRL, awsRR);
    display(k-31);
    x_store(k-31, :) = x;
    x1_store(k-31, :) = x1;
end

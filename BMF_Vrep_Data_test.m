clc;
close all;

% CONSTANTS
% filename = 'dataForMatlab_Moving.csv';
% filename = 'dataForMatlab_static.csv'
% filename = 'dataForMatlab_constAccel.csv';
filename = 'dataForMatlab_constVel.csv';

datafile = csvread(filename);
% Number of States
NN = 1000;
n_states=13;
if strcmp(filename,'dataForMatlab_Moving.csv') || strcmp(filename,'dataForMatlab_static.csv')
    n_meas_states = 7;
else
    n_meas_states = 8;
end
% Init state

s_state=[...
        datafile(32,2);...    %         datafile(32,2);...         % 1. x
        datafile(32,3);...    %         datafile(32,3);...         % 2. y
        datafile(32,11);...   %         datafile(32,11);...        % 3. vx
        datafile(32,12);...   %         datafile(32,12);...        % 4. vy
        datafile(32,14);...   %         datafile(32,14);...        % 5. ax
        datafile(32,15);...   %         datafile(32,15);...        % 6. ay
        datafile(32,5);...    %         datafile(32,5);...         % 7. yaw
        0;...                 % 8. yaw rate
        0;...                 % 9. yaw accel
        0;...                 % 10. alpha front
        0;...                 % 11. alpha rear
        0;...                 % 12. slip front
        0;...                 % 13. slip rear
        ];

    
% Q AND R
q_proc= 0.01;                                 %std of process 
r_meas= 0.01;                                 %std of measurement
% Q_proc=(q_proc^2)*eye(n_states);              % covariance of process
Q_proc = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.001 0.01 0.001 0.001 0.001 0.001 0.001]);
% R_meas=r_meas^2*eye(n_meas_states);

if strcmp(filename,'dataForMatlab_Moving.csv') || strcmp(filename,'dataForMatlab_static.csv')
    R_meas = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.001]);
else
    R_meas = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.001 0.001]);
end


    
x=s_state; %+q_proc*0.0001*randn(n_states,1);
P = diag([0.0 0.0 0.0 0.0 0.0 0.0 0.0 100 100 100 100 100 100]);
x_store = zeros(NN,n_states);
x1_store = zeros(NN,n_states);%estmate
rankObs_store = zeros(NN,n_states);
meas_store = zeros(NN, n_meas_states);             %act

sk = 110;
ek = NN;
if strcmp(filename,'dataForMatlab_Moving.csv') || strcmp(filename,'dataForMatlab_static.csv')
    dt_in = 0.13;
else
    dt_in = datafile(sk+1,1)-datafile(sk,1);
end




for k=sk: ek
    if strcmp(filename,'dataForMatlab_Moving.csv') || strcmp(filename,'dataForMatlab_static.csv')
        meas = [datafile(k, 2); datafile(k,3); datafile(k,11); datafile(k,12); datafile(k,14); datafile(k,15); datafile(k,5)];
    else
        meas = [datafile(k, 2); datafile(k,3); datafile(k,11); datafile(k,12); 0; 0; datafile(k,5);0];
    end
    delta = datafile(k,17);
    awsFL = datafile(k,18);
    awsFR = datafile(k,19);
    awsRL = datafile(k,20);
    awsRR = datafile(k,21);
    x_old = x;
    x_new = zeros(13);
    [x, P, x1] = BMF_ekf(x, P, Q_proc, R_meas, meas, delta, awsFL, awsFR, awsRL, awsRR, dt_in);
    x_new = x;
    if strcmp(filename,'dataForMatlab_Moving.csv') || strcmp(filename,'dataForMatlab_static.csv')
        dt_in = 0.13;
    else
        dt_in = datafile(k+1,1) - datafile(k,1);
    end
    
    display(k-sk+1);
    x_store(k-sk+1, :) = x;
    x1_store(k-sk+1, :) = x1;
    meas_store(k-sk+1, :) = meas;
%     rankObs_store(k-31, :) = rankObs;

end

figure;
plot(1:ek-sk, x_store(1:ek-sk,5),'b-',1:ek-sk, meas_store(1:ek-sk,5), 'r--');
figure;
plot(1:ek-sk, x_store(1:ek-sk,12),'b-');
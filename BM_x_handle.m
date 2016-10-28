function [xdash] = BM_x_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR)

dt = 0.130;
m = 1700; 	% Mass of vehicle in kgg
h = 0.84; 	% m
tf = 1.28;	% m
tr = tf;
%d = tf/2
a = 1.5;	% m
b = 1.5;	% m
g = 9.81;
Iz = 3825;
Ca_FL = 40000.0; %200520.0    # cornering stiffness / lateral stiffness
Ca_FR = 40000.0;
Ca_RL = 40000.0; %76426.0
Ca_RR = 40000.0;
Ck_FL = 150000.0;  %30.0        # braking stiffness / longitudinal stiffness
Ck_FR = 150000.0;  %30.0
Ck_RL = 150000.0;  %24.0
Ck_RR = 150000.0;  %24.0
wheelr = (0.4/2.0) + 1.1;

xdash = x;

xdash(0+1) = (x(0+1) + x(2+1)*dt + 0.5*x(4+1)*dt^2);  % x coord
xdash(1+1) = (x(1+1) + x(3+1)*dt + 0.5*x(5+1)*dt^2);  % y coord
xdash(2+1) = (x(2+1) + x(4+1) * dt);  % x vel
xdash(3+1) = (x(3+1) + x(5+1) * dt);  % y vel
xdash(4+1) = (x(4+1));  % x accel
xdash(5+1) = (x(5+1));  % y accel
xdash(6+1) = (x(6+1) + ((sqrt(x(2+1)^2 + x(3+1)^2))*tan(delta)*dt)/(a+b));  % yaw
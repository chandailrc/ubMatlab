function [xdash] = BMF_x_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in)

dt = dt_in;
m = 1300; 	% Mass of vehicle in kgg
h = 0.84; 	% m
tf = 1.436;	% m
tr = tf;
%d = tf/2
a = 1.0;	% m
b = 1.454;	% m
g = 9.81;
Iz = 1627;
Ca_FL = 30000.0;  % cornering stiffness / lateral stiffness
Ca_FR = 30000.0;
Ca_RL = 30000.0; 
Ca_RR = 30000.0;
Ck_FL = 50000.0;  % braking stiffness / longitudinal stiffness
Ck_FR = 50000.0;  
Ck_RL = 50000.0;  
Ck_RR = 50000.0;  
wheelr = 0.503; %(0.4/2.0) + 1.1;

xdash = x;

% Fxf = Ck_FL*x(11+1); Fyf = Ca_FL*x(9+1); Fxr = Ck_RL*x(12+1); Fyr = Ca_RL*x(10+1);

xdash(0+1) = (x(0+1) + x(2+1)*dt + 0.5*x(4+1)*dt^2);  % x coord
xdash(1+1) = (x(1+1) + x(3+1)*dt + 0.5*x(5+1)*dt^2);  % y coord
xdash(2+1) = (x(2+1) + x(4+1) * dt);  % x vel
xdash(3+1) = (x(3+1) + x(5+1) * dt);  % y vel
xdash(4+1) = ((1/m)*(Ck_FL*x(11+1)*cos(delta) - Ca_FL*x(9+1)*sin(delta) + Ck_RL*x(12+1))*cos(x(6+1))) - ((1/m)*(Ck_FL*x(11+1)*sin(delta) + Ca_FL*x(9+1)*cos(delta) + Ca_RL*x(10+1))*sin(x(6+1)));  % x accel
% display((1/m)*(Ck_FL*x(11+1)*cos(delta) - Ca_FL*x(9+1)*sin(delta) + Ck_RL*x(12+1))*cos(x(6+1)))
% display((1/m)*(Ck_FL*x(11+1)*sin(delta) + Ca_FL*x(9+1)*cos(delta) + Ca_RL*x(10+1))*sin(x(6+1)))
xdash(5+1) = ((1/m)*(Ck_FL*x(11+1)*sin(delta) + Ca_FL*x(9+1)*cos(delta) + Ca_RL*x(10+1))*cos(x(6+1))) + ((1/m)*(Ck_FL*x(11+1)*cos(delta) - Ca_FL*x(9+1)*sin(delta) + Ck_RL*x(12+1))*sin(x(6+1)));  % y accel
xdash(6+1) = (x(6+1) + ((sqrt(x(2+1)^2 + x(3+1)^2))*tan(delta)*dt)/(a+b));  % yaw
xdash(7+1) = (x(7+1) + x(8+1)*dt);  % yaw rate
xdash(8+1) = (1/Iz)*(a*(Ca_FL*x(9+1)*cos(delta) + Ck_FL*x(11+1)*sin(delta)) - b*(Ca_RL*x(10+1)));  % yaw accel
xdash(9+1) = (delta - (atan((x(3+1)*cos(x(6+1)) - x(2+1)*sin(x(6+1)) + a*x(7+1))/(x(3+1)*sin(x(6+1)) + x(2+1)*cos(x(6+1)))))); % alpha front
xdash(10+1) = -(atan((x(3+1)*cos(x(6+1)) - x(2+1)*sin(x(6+1)) - b*x(7+1))/(x(3+1)*sin(x(6+1)) + x(2+1)*cos(x(6+1)))));  % alpha rear
xdash(11+1) = -((((sqrt((x(3+1)*cos(x(6+1)) - x(2+1)*sin(x(6+1)) + a*x(7+1))^2 + (x(3+1)*sin(x(6+1)) + x(2+1)*cos(x(6+1)))^2))*cos(x(9+1)))/((angular_whlSpeed_FL+angular_whlSpeed_FR)*0.5*wheelr)) - 1);  % slip front
% display((sqrt((x(3+1)*cos(x(6+1)) - x(2+1)*sin(x(6+1)) + a*x(7+1))^2 + (x(3+1)*sin(x(6+1)) + x(2+1)*cos(x(6+1)))^2))*cos(x(9+1)));
% display((angular_whlSpeed_FL+angular_whlSpeed_FR)*0.5*wheelr);
xdash(12+1) = -((((sqrt((x(3+1)*cos(x(6+1)) - x(2+1)*sin(x(6+1)) - b*x(7+1))^2 + (x(3+1)*sin(x(6+1)) + x(2+1)*cos(x(6+1)))^2))*cos(x(10+1)))/((angular_whlSpeed_RL+angular_whlSpeed_RR)*0.5*wheelr)) - 1);  % slip rear


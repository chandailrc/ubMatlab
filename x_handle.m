function [xdash] = x_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR)

dt = 0.130;
m = 941.20; 	% Mass of vehicle in kgg
h = 0.84; 	% m
tf = 1.28;	% m
tr = tf;
%d = tf/2
a = 0.968;	% m
b = 1.392;	% m
g = 9.81;
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

xdash = x;

xdash(0+1) = (x(0+1) + x(2+1) * dt + 0.5 * x(4+1) * dt^2);  % x coord
xdash(1+1) = (x(1+1) + x(3+1) * dt + 0.5 * x(5+1) * dt^2);  % y coord
xdash(2+1) = (x(2+1) + x(4+1) * dt);  % x vel
xdash(3+1) = (x(3+1) + x(5+1) * dt);  % y vel
xdash(4+1) = (((1 / m) * (Ck_FL * x(24+1) * cos(delta) - Ca_FL * x(15+1) * sin(delta) + Ck_FR * x(25+1) * cos(delta) - Ca_FR * x(16+1) * sin(delta) + Ck_RL * x(26+1) + Ck_RR * x(27+1))) * cos(x(6+1)) - ((1 / m) * (Ck_FL * x(24+1) * sin(delta) + Ca_FL * x(15+1) * cos(delta) + Ck_FR * x(25+1) * sin(delta) + Ca_FR * x(16+1) * cos(delta) + Ca_RL * x(17+1) + Ca_RR * x(18+1))) * sin(x(6+1)));  % x accel
xdash(5+1) = (((1 / m) * (Ck_FL * x(24+1) * sin(delta) + Ca_FL * x(15+1) * cos(delta) + Ck_FR * x(25+1) * sin(delta) + Ca_FR * x(16+1) * cos(delta) + Ca_RL * x(17+1) + Ca_RR * x(18+1))) * cos(x(6+1)) + ((1 / m) * (Ck_FL * x(24+1) * cos(delta) - Ca_FL * x(15+1) * sin(delta) + Ck_FR * x(25+1) * cos(delta) - Ca_FR * x(16+1) * sin(delta) + Ck_RL * x(26+1) + Ck_RR * x(27+1))) * sin(x(6+1)));  % y accel
xdash(6+1) = (x(6+1) + x(9+1) * dt + 0.5 * x(12+1) * dt^2);  % yaw
xdash(7+1) = x(7+1);  % pitch
xdash(8+1) = x(8+1);  % roll
xdash(9+1) = x(9+1) + x(12+1) * dt;  % yaw rate
xdash(10+1) = x(10+1);  % pitch rate
xdash(11+1) = x(11+1);  % roll rate
xdash(12+1) = (1 / Iz) * (-tf * 0.5 * (Ck_FL * x(24+1) * cos(delta) - Ca_FL * x(15+1) * sin(delta)) + tf * 0.5 * (Ck_FR * x(25+1) * cos(delta) - Ca_FR * x(16+1) * sin(delta)) - tr * 0.5 * Ck_RL * x(26+1) + tr * 0.5 * Ck_RR * x(27+1) + a * (Ca_FL * x(15+1) * cos(delta) + Ck_FL * x(24+1) * sin(delta)) + a * (Ca_FR * x(16+1) * cos(delta) + Ck_FR * x(25+1) * sin(delta)) - b * Ca_RL * x(17+1) - b * Ca_RR *x(18+1));  % yaww accel
%display(-tf * 0.5 * (Ck_FL * x(24+1) * cos(delta) - Ca_FL * x(15+1) * sin(delta)));
%display(tf * 0.5 * (Ck_FR * x(25+1) * cos(delta) - Ca_FR * x(16+1) * sin(delta)));
%display(- tr * 0.5 * Ck_RL * x(26+1) + tr * 0.5 * Ck_RR * x(27+1));
%display(a * (Ca_FL * x(15+1) * cos(delta) + Ck_FL * x(24+1) * sin(delta)));
%display(a * (Ca_FR * x(16+1) * cos(delta) + Ck_FR * x(25+1) * sin(delta)));
%display(- b * Ca_RL * x(17+1) - b * Ca_RR *x(18+1));
%display((1 / Iz) * (-tf * 0.5 * (Ck_FL * x(24+1) * cos(delta) - Ca_FL * x(15+1) * sin(delta)) + tf * 0.5 * (Ck_FR * x(25+1) * cos(delta) - Ca_FR * x(16+1) * sin(delta)) - tr * 0.5 * Ck_RL * x(26+1) + tr * 0.5 * Ck_RR * x(27+1) + a * (Ca_FL * x(15+1) * cos(delta) + Ck_FL * x(24+1) * sin(delta)) + a * (Ca_FR * x(16+1) * cos(delta) + Ck_FR * x(25+1) * sin(delta)) - b * Ca_RL * x(17+1) - b * Ca_RR *x(18+1)));
xdash(13+1) = x(13+1);  % pitch accel
xdash(14+1) = x(14+1);  % roll accel
xdash(15+1) = delta - atan((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) + a * x(9+1)) / (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) - tf * 0.5 * x(9+1)));  % alpha slip angle FL + ve
xdash(16+1) = delta - atan((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) + a * x(9+1)) / (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) + tf * 0.5 * x(9+1)));  % alpha slip angle FR + ve
xdash(17+1) = atan((-(x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1))) + b * x(9+1)) / (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) - tr * 0.5 * x(9+1)));  % alpha slip angle RL -ve
xdash(18+1) = atan((-(x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1))) + b * x(9+1)) / (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) + tr * 0.5 * x(9+1)));  % alpha slip angle RR -ve
xdash(19+1) = atan((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1))) / (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1))));  % Beta slip angle
xdash(20+1) = sqrt((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) + a * x(9+1))^2 + (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) - tf * 0.5 * x(9+1))^2) * cos(x(15+1));  % Actual wheel velocity FL
xdash(21+1) = sqrt((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) + a * x(9+1))^2 + (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) + tf * 0.5 * x(9+1))^2) * cos(x(16+1));  % Actual wheel velocity FR
xdash(22+1) = sqrt((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) - b * x(9+1))^2 + (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) - tr * 0.5 * x(9+1))^2) * cos(x(17+1));  % Actual wheel velocity RL
xdash(23+1) = sqrt((x(3+1) * cos(x(6+1)) - x(2+1) * sin(x(6+1)) - b * x(9+1))^2 + (x(3+1) * sin(x(6+1)) + x(2+1) * cos(x(6+1)) + tr * 0.5 * x(9+1))^2) * cos(x(18+1));  % Actual wheel velocity RR
xdash(24+1) = 1 - (x(20+1)/(angular_whlSpeed_FL * wheelr));	% slip ratio FL
xdash(25+1) = 1 - (x(21+1)/(angular_whlSpeed_FR * wheelr));  % slip ratio FR
xdash(26+1) = 1 - (x(22+1)/(angular_whlSpeed_RL * wheelr));  % slip ratio RL
xdash(27+1) = 1 - (x(23+1)/(angular_whlSpeed_RR * wheelr));  % slip ratio RR
xdash(28+1) = 1.0; %(0.5 * m * g - m * x(5+1) * h / tf) * b / (a + b) - 0.5 * m * x(4+1) * h / (a + b);  % Fz FL
xdash(29+1) = 1.0; %(0.5 * m * g + m * x(5+1) * h / tf) * b / (a + b) - 0.5 * m * x(4+1) * h / (a + b);  % Fz FR
xdash(30+1) = 1.0; %(0.5 * m * g - m * x(5+1) * h / tr) * a / (a + b) + 0.5 * m * x(4+1) * h / (a + b);  % Fz RL
xdash(31+1) = 1.0; %(0.5 * m * g + m * x(5+1) * h / tr) * a / (a + b) + 0.5 * m * x(4+1) * h / (a + b);  % Fz RR


function [ydash] = BM_y_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR)
ydash(1) = x(0+1);
ydash(2) = x(1+1);
ydash(3) = x(2+1);
ydash(4) = x(3+1);
ydash(5) = x(4+1);
ydash(6) = x(5+1);
ydash(7) = x(6+1);
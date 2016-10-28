function [x_out , P_out] = ekf(x_in, P_in, Q, R, meas, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR)
x1 = x_handle(x_in, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR);
A = x_Jacob(x_in, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR);
P=A*P_in*A'+Q;
z1 = (y_handle(x1))';
H = y_Jacob();
P12=P*H';

% eig_A = eig(H*P12+R);
% flag = 0;
% for i = 1:rank(H*P12+R)
% 	if eig_A(i) <= 0 
% 	flag = 1;
% 	end
% end


R=chol(H*P12+R);            %Cholesky factorizationde
U=P12/R;                    %K=U/R'; Faster because of back substitution
x_out=x1+U*(R'\(meas-z1));         %Back substitution to get state update
P_out=P-U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12.

%K=P12/(H*P12+R);
%x_out=x1+K*(meas-z1);
%P_out=P-K*P12';

function [x,P,x1]=BMF_ekf(x,P,Q,R,z, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in)

[x1,A]=jaccsd(1,x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in);    %nonlinear update and linearization at current state
P=A*P*A'+Q;                 %partial update
[z1,H]=jaccsd(2,x1, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in);    %nonlinear measurement and linearization
P12=P*H';                   %cross covariance
retRankObs = rank(obsv(A,H)); 
display(retRankObs);
% display(obsv(A,H));

% K=P12*inv(H*P12+R);       %Kalman filter gain
% x=x1+K*(z-z1);            %state estimate
% P=P-K*P12';               %state covariance matrixS
R=chol(H*P12+R);            %Cholesky factorization
U=P12/R;                    %K=U/R'; Faster because of back substitution
x=x1+U*(R'\(z-z1'));         %Back substitution to get state update
P=P-U*U';                   %Covariance update, U*U'=P12/R/R'*P12'=K*P12.

function [z,A]=jaccsd(chkflag,x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
if chkflag == 1
    z=BMF_x_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in);
    n=numel(x);
    m=numel(z);
    A=zeros(m,n);
    h=n*eps;
    for k=1:n
        x1=x;
        x1(k)=x1(k)+h*i;
        A(:,k)=imag(BMF_x_handle(x1, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in))/h;
    end
else
    z=BMF_y_handle(x, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in);
    n=numel(x);
    m=numel(z);
    A=zeros(m,n);
    h=n*eps;
    for k=1:n
        x1=x;
        x1(k)=x1(k)+h*i;
        A(:,k)=imag(BMF_y_handle(x1, delta, angular_whlSpeed_FL, angular_whlSpeed_FR, angular_whlSpeed_RL, angular_whlSpeed_RR, dt_in))/h;
    end
end
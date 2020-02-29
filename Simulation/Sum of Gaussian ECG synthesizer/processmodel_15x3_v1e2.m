function y=processmodel_15x3_v1e2(x,delta,beta,theta)

% y: the next samples of state variables
% x: the state variables vector at the current sample
% delta: sampling period
% beta: linear coefficient of angular velocity 
% theta: 15_dim parameters vector, containing 
%          a_P, a_Q, a_R, a_S, a_T,
%          b_P, b_Q, b_R, b_S, b_T,
%          c_P, c_Q, c_R, c_S, c_T,
% w: Process noise 
%          eta=process noise added to quassian modelling functions (z)
%          u= process noise for angular velocity





%%% p wave
dc(1)=rem((x(1)-theta(11)+pi),2*pi)-pi; % distance from the center
ep(1)=exp(-(dc(1)^2)./(2*theta(6)^2)); % Gaussian function (differential form) 
zp=(delta*theta(1)*x(3))*dc(1)*(ep(1))./(theta(6)^2); % scaling the Gaussian (differential form)


%%% Q wave
dc(2)=rem((x(1)-theta(12)+pi),2*pi)-pi;
ep(2)=exp(-(dc(2)^2)./(2*theta(7)^2));
zq=(delta*theta(2)*x(3))*dc(2)*(ep(2))./(theta(7)^2);

%%% R wave
dc(3)=rem((x(1)-theta(13)+pi),2*pi)-pi;
ep(3)=exp(-(dc(3)^2)./(2*theta(8)^2));
zr=(delta*theta(3)*x(3))*dc(3)*(ep(3))./(theta(8)^2);

%%% S wave
dc(4)=rem((x(1)-theta(14)+pi),2*pi)-pi;
ep(4)=exp(-(dc(4)^2)./(2*theta(9)^2));
zs=(delta*theta(4)*x(3))*dc(4)*(ep(4))./(theta(9)^2);

%%% T wave
dc(5)=rem((x(1)-theta(15)+pi),2*pi)-pi;
ep(5)=exp(-(dc(5)^2)./(2*theta(10)^2));
zt=(delta*theta(5)*x(3))*dc(5)*(ep(5))./(theta(10)^2);


y=-(zp+zq+zr+zs+zt)+0.99*x(2); % updating the synthetic ECG
end










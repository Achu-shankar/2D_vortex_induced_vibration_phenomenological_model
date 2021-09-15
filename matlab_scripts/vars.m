% Vr = 10.9;
CDo = 0.2;
CLo = 0.3;
St  = 0.2;
rho = 1000;
D = 0.114;
CM = 1;   
ms = 0.1;  % cylinder mass
gamma = 0.8;

zeta_x = 0.002;
zeta_y = 0.002;

freq_rat = 1;

mf = pi*rho*D^2*CM/4; % fluid added mass
omega = St*Vr;

% mu  = (ms+mf)/(rho*D^2);
% m_str = 4*mu/pi-CM;
% if t<0.002
%     mu
%     m_str
% end
m_str = 5.4;
mu   = (m_str + CM)*pi/4;
MD  = CDo/(16*pi^2*St^2*mu);
ML  = CLo/(16*pi^2*St^2*mu);



lambda_x = 2*zeta_x*freq_rat + gamma*omega/mu;
lambda_y = 2*zeta_y + gamma*omega/mu;    

alpha_x = 0.7;
beta_x = 0.7;
eps_x = 0.3;
LAMBDA_x = 15;

alpha_y = 0.7;
beta_y = 0.7;
eps_y = 0.00234*exp(0.228*m_str);  % for freq_rat = 1
% eps_y = 0.0045*exp(0.228*m_str);  % for freq_rat = 1
LAMBDA_Y =15;


%% 2D VIV 
% close all
clearvars
% X  = 1-x,2-x_d,3-p,4-p_d,5-y,6-y_d,7-q,8-q_d


%% Variation of max amplitude with reduced velocity 
% params
X_init = [0;0;2;0;0;0;2;0];
time = [0,100];
% Vr = 1:0.5:2.5;
Vr = 10.9;
delt = 0.1;
Xmax = zeros(size(Vr));
Ymax = zeros(size(Vr));
Xrms = zeros(size(Vr));
Yrms = zeros(size(Vr));
Ntl   = 100/delt ;
Ntf  = Ntl*0.7;
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
foy = zeros(size(Vr));
fox = zeros(size(Vr));

for i = 1:length(Vr)
    [t,sol] = RK4(@(t,X)res(t,X,Vr(i)),time,delt,X_init);
%      [t,sol] = ode45(@(t,X)res(t,X,Vr(i)),time,X_init,opts);
%      sol = sol';
%      Ntl  = length(t)-2 ;
%      Ntf  = int64(Ntl*0.7);
     sol(1,:) = sol(1,:) - (max(sol(1,Ntf:Ntl))+ min(sol(1,Ntf:Ntl)))/2 ;
     sol(5,:) = sol(5,:) - (max(sol(5,Ntf:Ntl)) + min(sol(5,Ntf:Ntl)))/2 ;
     
     Fs = Ntl-Ntf;
     n = 2^nextpow2(Fs+1);
     Y = fft(sol(1,Ntf:Ntl),n);
     f = Fs*(0:(n/2))/(n*time(2));
     P = abs(Y/n);
     [val,ind] = max(P);
     fox(i) = 2*f(ind);
     n = 2^nextpow2(Fs+1);
     Y = fft(sol(5,Ntf:Ntl),n);
     f = Fs*(0:(n/2))/(n*time(2));
     P = abs(Y/n);
     [val,ind] = max(P);
     foy(i) = 2*f(ind);
          
     Xmax(i) = max(sol(1,Ntf:Ntl));
     Ymax(i) = max(sol(5,Ntf:Ntl));
     Xrms(i) = rms(sol(1,Ntf:Ntl));
     Yrms(i) = rms(sol(5,Ntf:Ntl));
     i
end
foy = foy/0.3120;
fox = fox/0.316;

%% 
figure
subplot(2,1,1)
plot(t,sol(1,:))
ylim([-1,1])
xlim([50,100])
grid on
title('X inline')
% 
subplot(2,1,2)
plot(t,sol(5,:))
ylim([-2,2])
xlim([50,100])
grid on
title('Y cross flow')

figure
plot(sol(1,Ntf:Ntl),sol(5,Ntf:Ntl))
figure
plot(Vr,Xmax,'r','linewidth',2)
hold on
plot(Vr,Xrms,'b','linewidth',2)
title('X')
ylim([0,0.7])
figure
plot(Vr,Ymax,'r','linewidth',2)
hold on
plot(Vr,Yrms,'b','linewidth',2)
title('Y')
ylim([0,2])
% subroutine 

function dX=res(t,X,Vr)
% params
% consider theta is clockwise

vars
dx   = X(2);
dx_d = - lambda_x*X(2) - freq_rat^2*(X(1)+alpha_x*X(1)^3 ...
       + beta_x*X(1)*X(5)^2) + MD*omega^2*X(3) - 2*pi*ML*omega^2*X(7)*X(6)/Vr;
   
dp   = X(4);
dp_d = -2*eps_x*omega*(X(3)^2-1)*X(4)  - 4*omega^2*X(3) + LAMBDA_x*dx_d;

dy   = X(6);
dy_d = - lambda_y*X(6) - X(5) - alpha_y*X(5)^3 - beta_y*X(5)*X(1)^2 ...
       + ML*omega^2*X(7) + 2*pi*MD*omega^2*X(3)*X(6)/Vr;
   
dq   = X(8);
dq_d = -eps_y*omega*(X(7)^2-1)*X(8) - omega^2*X(7) + LAMBDA_Y*dy_d;

dX = [dx;dx_d;dp;dp_d;dy;dy_d;dq;dq_d];

end


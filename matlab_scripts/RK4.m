%% Runge Kutta 4th order ode solver
function [t,sol] = RK4(func,time,delt,init)
    t = time(1):delt:time(2);
    X = zeros(length(init),length(t));
    X(:,1) = init;
    for i=1:length(t)-1
        k1 = delt*func(t(i),X(:,i));
        k2 = delt*func(t(i)+delt/2,X(:,i)+k1/2);
        k3 = delt*func(t(i)+delt/2,X(:,i)+k2/2);
        k4 = delt*func(t(i)+delt,X(:,i)+k3);
        X(:,i+1) = X(:,i) + (k1+2*k2+2*k3+k4)/6;
    end
    sol = X;
end
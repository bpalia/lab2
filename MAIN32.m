clearvars
close all

load tork.dat

cf = 40;

u = tork(:,2); %input
y = tork(:,1); %output
u = u - mean(u);
y = y - mean(u);
z = iddata(y,u)

fnum = 0;
fnum = fnum+1;
figure(fnum)
plot(z(1:300))

%% prewithening u
fnum = fnum+1;
figure(fnum)
subplot(211)
acf(u, cf, 0.05, true, 0, 0);
title('ACF for u')
subplot(212)
pacf(u, cf, 0.05, true, 0);
title('PACF for u')
fnum = fnum + 1;
figure(fnum)
normplot(u)

data_u = iddata(u);
ar_model_u = arx(data_u, 1);
present(ar_model_u)
res_u = resid(ar_model_u, data_u);

fnum = fnum+1;
figure(fnum)
subplot(211)
acf(res_u.y, cf, 0.05, true, 0, 0);
title('ACF for resid AR1 for u')
subplot(212)
pacf(res_u.y, cf, 0.05, true, 0);
title('PACF for resid AR1 for u')
fnum = fnum + 1;
figure(fnum)
normplot(res_u.y)

upw = res_u.y;
data_y = iddata(y);
res_y = resid(ar_model_u, data_y);
ypw = res_y.y;
%%
fnum = fnum +1;
figure(fnum)
stem(-cf:cf,crosscorr(upw,ypw,cf));
title('Cross correlation function u_{pw} and y_{pw}'); xlabel('Lag');
hold on
plot(-cf:cf,2/sqrt(length(upw))*ones(1,2*cf+1),'--') 
plot(-cf:cf,-2/sqrt(length(upw))*ones(1,2*cf+1),'--')
hold off
% d = 3, r = 2, s = 1
%% estimating B and A2
A2 = [1 0 0];
B = [0 0 0 0 0];
Mi = idpoly(1,B,[],[],A2);
Mi.Structure.b.Free = [zeros(1, 3), ones(1,2)];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); present(Mba2);
vhat = resid(Mba2,zpw);

fnum = fnum +1;
figure(fnum)
stem(-cf:cf,crosscorr(upw,vhat.y,cf));
title('Cross correlation function u_{pw} and v'); xlabel('Lag');
hold on
plot(-cf:cf,2/sqrt(length(upw))*ones(1,2*cf+1),'--') 
plot(-cf:cf,-2/sqrt(length(upw))*ones(1,2*cf+1),'--')
hold off

%% estimating C1 and A1
h = filter(Mba2.b, Mba2.f,u);
x = y(51:end) - h(51:end);

fnum = fnum +1;
figure(fnum)
stem(-cf:cf,crosscorr(u(51:end),x,cf));
title('Cross correlation function u and x'); xlabel('Lag');
hold on
plot(-cf:cf,2/sqrt(length(x))*ones(1,2*cf+1),'--') 
plot(-cf:cf,-2/sqrt(length(x))*ones(1,2*cf+1),'--')
hold off

fnum = fnum+1;
figure(fnum)
subplot(211)
acf(x, cf, 0.05, true, 0, 0);
title('ACF for x')
subplot(212)
pacf(x, cf, 0.05, true, 0);
title('PACF for x')
fnum = fnum + 1;
figure(fnum)
normplot(x)

% Estimating AR1
data_x = iddata(x);
ar1_model_x = arx(data_x, 1);
present(ar1_model_x)
res_x = resid(ar1_model_x,data_x);
fnum = fnum +1;
figure(fnum)
subplot(211)
acf(res_x.y, cf, 0.05, true, 0, 0);
title('ACF for resid of x AR1')
subplot(212)
pacf(res_x.y, cf, 0.05, true, 0);
title('PACF for resid of x AR1')
fnum = fnum + 1;
figure(fnum)
normplot(res_x.y)

% We have AR model C = 1;
%% All the polynomials
A1 = [1 0];
A2 = [1 0 0];
B = [0 0 1 0 0];
C = [1];
Mi = idpoly(1,B,C,A1,A2);
%Mi.Structure.b.Free = [zeros(1, 3), 1, 1];
z = iddata(y,u);
MboxJ = pem(z,Mi); present(MboxJ);
ehat = resid(MboxJ,z);

fnum = fnum +1;
figure(fnum)
stem(-cf:cf,crosscorr(u,ehat.y,cf));
title('Cross correlation function u and e'); xlabel('Lag');
hold on
plot(-cf:cf,2/sqrt(length(u))*ones(1,2*cf+1),'--') 
plot(-cf:cf,-2/sqrt(length(u))*ones(1,2*cf+1),'--')
hold off

fnum = fnum +1;
figure(fnum)
subplot(211)
acf(ehat.y, cf, 0.05, true, 0, 0);
title('ACF for resid e')
subplot(212)
pacf(ehat.y, cf, 0.05, true, 0);
title('PACF for resid e')
fnum = fnum + 1;
figure(fnum)
normplot(ehat.y)
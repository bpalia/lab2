clearvars
close all

fnum = 0;
cutoff = 50;

% n = 500;
% A1 = [1 -.65]; A2 = [ 1 .90 .78];
% C = 1; B = [0 0 0 0 .4];
% e = sqrt(1.5) * randn(n + 100, 1);
% w = sqrt(2) * randn(n + 200, 1 );
% A3 = [1 .5]; C3 = [1 -.3 .2];
% u = filter(C3,A3,w); u = u(101:end);
% y = filter(C,A1,e) + filter(B,A2,u);
% u = u(101:end); y = y(101:end);
% clear A1, A2, C, B, e, w, A3, C3
load data_lab2.mat

fnum = fnum+1;
figure(fnum)
subplot(211)
acf(u, cutoff, 0.05, true, 0, 0);
title('ACF for u')
subplot(212)
pacf(u, cutoff, 0.05, true, 0);
title('PACF for u')
fnum = fnum + 1;
figure(fnum)
normplot(u)

data_u = iddata(u);
%% AR 1
ar_model_u = arx(data_u, 1);
present(ar_model_u)
res1 = filter(ar_model_u.a, 1, u);
res1 = res1(51:end);

fnum = fnum+1;
figure(fnum)
subplot(211)
acf(res1, cutoff, 0.05, true, 0, 0);
title('ACF for resid AR1')
subplot(212)
pacf(res1, cutoff, 0.05, true, 0);
title('PACF for resid AR1')
fnum = fnum + 1;
figure(fnum)
normplot(res1)
fnum = fnum+1;
figure(fnum)
whitenessTest(res1, 0.05)

%% ARMA11
arma_model11_u = armax(data_u, [1 1]);
present(arma_model11_u)

res2 = filter(arma_model11_u.a, arma_model11_u.c, u);
res2 = res2(51:end);
fnum = fnum+1;
figure(fnum)
subplot(211)
acf(res2, cutoff, 0.05, true, 0, 0);
title('ACF for resid ARMA11')
subplot(212)
pacf(res2, cutoff, 0.05, true, 0);
title('PACF for resid ARMA11')
fnum = fnum + 1;
figure(fnum)
normplot(res2)
fnum = fnum+1;
figure(fnum)
whitenessTest(res2, 0.05)

%% ARMA12
arma_model12_u = armax(data_u, [1 2]);
present(arma_model12_u)

res3 = filter(arma_model12_u.a, arma_model12_u.c, u);
res3 = res3(51:end);
fnum = fnum+1;
figure(fnum)
subplot(211)
acf(res3, cutoff, 0.05, true, 0, 0);
title('ACF for resid ARMA12')
subplot(212)
pacf(res3, cutoff, 0.05, true, 0);
title('PACF for resid ARMA12')
fnum = fnum + 1;
figure(fnum)
normplot(res3)
fnum = fnum+1;
figure(fnum)
whitenessTest(res3, 0.05)

%% Whitening
upw = res3;
fnum = fnum +1;
figure(fnum)
subplot(211)
acf(y, cutoff, 0.05, true, 0, 0);
title('ACF for y')
subplot(212)
pacf(y, cutoff, 0.05, true, 0);
title('PACF for y')
fnum = fnum + 1;
figure(fnum)
normplot(y)

ypw = filter(arma_model12_u.a, arma_model12_u.c, y);
ypw = ypw(51:end);

%% 
fnum = fnum +1;
figure(fnum)
stem(-cutoff:cutoff,crosscorr(upw,ypw,cutoff));
title('Cross correlation function u_{pw} and y_{pw}'); xlabel('Lag');
hold on
plot(-cutoff:cutoff,2/sqrt(length(upw))*ones(1,2*cutoff+1),'--') % 700 samples as for w
plot(-cutoff:cutoff,-2/sqrt(length(upw))*ones(1,2*cutoff+1),'--')
hold off

%% incorrect d = 4; r = 0; s = 0;
A2 = [1];
B = [0 0 0 0 0];
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [zeros(1, 4), 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); present(Mba2);
vhat = resid(Mba2,zpw);

fnum = fnum +1;
figure(fnum)
stem(-cutoff:cutoff,crosscorr(upw,vhat.y,cutoff));
title('Cross correlation function u_{pw} and v'); xlabel('Lag');
hold on
plot(-cutoff:cutoff,2/sqrt(length(upw))*ones(1,2*cutoff+1),'--') % 700 samples as for w
plot(-cutoff:cutoff,-2/sqrt(length(upw))*ones(1,2*cutoff+1),'--')
hold off

%% correct d = 4; r = 2; s = 0;
A2 = [1 0 0];
B = [0 0 0 0 0];
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [zeros(1, 4), 1];
zpw = iddata(ypw,upw);
Mba2 = pem(zpw,Mi); present(Mba2);
vhat = resid(Mba2,zpw);

fnum = fnum +1;
figure(fnum)
stem(-cutoff:cutoff,crosscorr(upw,vhat.y,cutoff));
title('Cross correlation function u_{pw} and v correct'); xlabel('Lag');
hold on
plot(-cutoff:cutoff,2/sqrt(length(upw))*ones(1,2*cutoff+1),'--') % 700 samples as for w
plot(-cutoff:cutoff,-2/sqrt(length(upw))*ones(1,2*cutoff+1),'--')
hold off

fnum = fnum +1;
figure(fnum)
subplot(211)
acf(vhat.y, cutoff, 0.05, true, 0, 0);
title('ACF for v')
subplot(212)
pacf(vhat.y, cutoff, 0.05, true, 0);
title('PACF for v')
fnum = fnum + 1;
figure(fnum)
normplot(vhat.y)

close all
clearvars
fnum = 0;
cf = 50;

load svedala
y = svedala;

fnum = fnum +1;
figure(fnum)
subplot(211)
acf(y, cf, 0.05, true, 0, 0);
title('ACF for data')
subplot(212)
pacf(y, cf, 0.05, true, 0);
title('PACF for data')
fnum = fnum + 1;
figure(fnum)
normplot(y)

% Deseasonalizing
A24 = [1, zeros(1,23), -1];
y_s = filter(A24, 1, y);
y_s(1:24) = [];

fnum = fnum +1;
figure(fnum)
subplot(211)
acf(y_s, cf, 0.05, true, 0, 0);
title('ACF for data S=24')
subplot(212)
pacf(y_s, cf, 0.05, true, 0);
title('PACF for data S=24')
fnum = fnum + 1;
figure(fnum)
normplot(y_s)

data = iddata(y_s);
model_init = idpoly([1, 0, 0], [], [1, zeros(1, 24)]);
model_init.Structure.c.Free = [zeros(1,24), 1];
model_armax = pem(data, model_init);

present(model_armax)
y_f = filter(model_armax.a, model_armax.c, y_s);

fnum = fnum +1;
figure(fnum)
subplot(211)
acf(y_f, cf, 0.05, true, 0, 0);
title('ACF for data ARMA(2,24) S=24')
subplot(212)
pacf(y_f, cf, 0.05, true, 0);
title('PACF for data ARMA(2,24) S=24')
fnum = fnum + 1;
figure(fnum)
normplot(y_f)
%%
A1 = model_armax.a;
C = model_armax.c;

A = conv(A1,A24);
%%
k = 3;
[F,G] = func_poldiv(A,C,k);
yhat3 = filter(G,C,y);

fnum = fnum+1;
figure(fnum)
plot(y)
title('3-step prediction')
hold on
plot(yhat3)
hold off

err3 = y(k+1:end)-yhat3(k+1:end);
err3_var = var(err3);

%%
k = 26;
[F,G] = func_poldiv(A,C,k);
yhat26 = filter(G,C,y);

fnum = fnum+1;
figure(fnum)
plot(y)
title('26-step prediction')
hold on
plot(yhat26)
hold off

err26 = y(k+1:end)-yhat26(k+1:end);
err26_var = var(err26);
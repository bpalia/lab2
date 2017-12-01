close all
clearvars
clc

fnum = 0;
load sturup
load svedala
y = svedala;
x = sturup;

A = [1 -1.49 0.57];
B = [0 0 0 0.28 -0.26];
C = [1];

%% k = 3;
k = 3;
[F,G] = func_poldiv(A,C,k);
BF = conv(B,F);
[Fhat,Ghat] = func_poldiv(C,BF,k);

yhat3 = filter(Ghat,C,x) + filter(G,C,y) + filter(Fhat,1,x);

fnum = fnum+1;
figure(fnum)
plot(y)
title('3-step prediction')
hold on
plot(yhat3)
hold off

err3 = y(k+1:end)-yhat3(k+1:end);
var_err3 = var(err3);

fnum = fnum +1;
figure(fnum)
plot(err3);
title('Prediction errors k = 3');

%% k = 26;
k = 26;
[F,G] = func_poldiv(A,C,k);
BF = conv(B,F);
[Fhat,Ghat] = func_poldiv(C,BF,k);

yhat26 = filter(Ghat,C,x) + filter(G,C,y) + filter(Fhat,1,x);

fnum = fnum+1;
figure(fnum)
plot(y)
title('26-step prediction')
hold on
plot(yhat26)
hold off

err26 = y(k+1:end)-yhat26(k+1:end);
var_err26 = var(err26);

fnum = fnum +1;
figure(fnum)
plot(err26);
title('Prediction errors k = 26');



clearvars
close all

load svedala
y = svedala;

A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];

[CS,AS] = equalLength(C,A);
k = 1;
[Fk, Gk] = deconv(conv([1, zeros(1,k-1)],CS),AS);

yhat_1 = filter(Gk,C,y);

fnum = 0;
fnum = fnum+1;
figure(fnum)
plot(y)
title('1-step prediction')
hold on
plot(yhat_1)
hold off

ehat_1 = y(2:end)-yhat_1(2:end);
var_n = var(ehat_1);
%% k = 3
k = 3;
[Fk, Gk] = deconv(conv([1, zeros(1,k-1)],CS),AS);
yhat_3 = filter(Gk,C,y);

fnum = fnum+1;
figure(fnum)
plot(y)
title('3-step prediction')
hold on
plot(yhat_3)
hold off

err_3 = y(k+1:end)-yhat_3(k+1:end);
err_3_m = mean(err_3);
err_3_var = sum(Fk.^2)*var_n;
CI_3 = 2*sqrt(var_n*sum(Fk.^2));

fnum = fnum +1;
figure(fnum)
plot(err_3);
title('Prediction errors k = 3');
hold on
plot(CI_3*ones(1,length(err_3)),'--') 
plot(-CI_3*ones(1,length(err_3)),'--')
hold off

numout_3 = sum(abs(err_3)>CI_3)/length(err_3)*100;

fnum = fnum+1;
figure(fnum)
R = covf(err_3,40);
stem(R) title ('Covariance of residuals for 3-step prediction')
%% k = 26
k = 26;
[Fk, Gk] = deconv(conv([1, zeros(1,k-1)],CS),AS);
yhat_26 = filter(Gk,C,y);

fnum = fnum+1;
figure(fnum)
plot(y)
title('26-step prediction')
hold on
plot(yhat_26)
hold off

err_26 = y(k+1:end)-yhat_26(k+1:end);
err_26_m = mean(err_26);
err_26_var = sum(Fk.^2)*var_n;
CI_26 = 2*sqrt(var_n*sum(Fk.^2));

fnum = fnum +1;
figure(fnum)
plot(err_26);
title('Prediction errors k = 26');
hold on
plot(CI_26*ones(1,length(err_26)),'--') 
plot(-CI_26*ones(1,length(err_26)),'--')
hold off

numout_26 = sum(abs(err_26)>CI_26)/length(err_26)*100;

fnum = fnum+1;
figure(fnum)
R = covf(err_26,40);
stem(R); title ('Covariance of residuals for 26-step prediction')
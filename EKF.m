clc
clear

% Init sistem
delT=0.001;
tcons=0.5;
varW=0.4;
xinit = [ones(1,7)]';
I = eye(7);

% Sistem param
A = [exp(-delT/tcons) zeros(1,6);
    0 exp(-delT/tcons) zeros(1,5);
    0 0 exp(-delT/tcons) zeros(1,4);
    (-(xinit(5)*delT)/2) (-(xinit(6)*delT)/2) (-(xinit(7)*delT)/2) 1 (-(xinit(1)*delT)/2) (-(xinit(2)*delT)/2) (-(xinit(3)*delT)/2);
    ((xinit(4)*delT)/2) (-(xinit(7)*delT)/2) ((xinit(6)*delT)/2) ((xinit(1)*delT)/2) 1 ((xinit(3)*delT)/2) (-(xinit(2)*delT)/2);
    ((xinit(7)*delT)/2) ((xinit(4)*delT)/2) (-(xinit(5)*delT)/2) ((xinit(2)*delT)/2) (-(xinit(3)*delT)/2) 1 ((xinit(1)*delT)/2);
    (-(xinit(6)*delT)/2) ((xinit(5)*delT)/2) ((xinit(4)*delT)/2) ((xinit(3)*delT)/2) ((xinit(2)*delT)/2) (-(xinit(1)*delT)/2) 1];
H = I;
Q = [eye(3)*(varW/(2*tcons)*1-exp(-delT/tcons)) zeros(3,4);zeros(4,7)];
R = [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001].*eye(7);

%% get data awal reference dengan QUEST algorithm
xk = xinit;
Adat = A;
for i = 1:200
    xk1 = Adat*xk + [(eye(3)*varW) zeros(3,4);zeros(4,7)]*randn(7,1);
    zk = H*xk + [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001]*eye(7)*randn(7,1);
    Adat = [exp(-delT/tcons) zeros(1,6);
        0 exp(-delT/tcons) zeros(1,5);
        0 0 exp(-delT/tcons) zeros(1,4);
        (-(xk1(5)*delT)/2) (-(xk1(6)*delT)/2) (-(xk1(7)*delT)/2) 1 (-(xk1(1)*delT)/2) (-(xk1(2)*delT)/2) (-(xk1(3)*delT)/2);
        ((xk1(4)*delT)/2) (-(xk1(7)*delT)/2) ((xk1(6)*delT)/2) ((xk1(1)*delT)/2) 1 ((xk1(3)*delT)/2) (-(xk1(2)*delT)/2);
        ((xk1(7)*delT)/2) ((xk1(4)*delT)/2) (-(xk1(5)*delT)/2) ((xk1(2)*delT)/2) (-(xk1(3)*delT)/2) 1 ((xk1(1)*delT)/2);
        (-(xk1(6)*delT)/2) ((xk1(5)*delT)/2) ((xk1(4)*delT)/2) ((xk1(3)*delT)/2) ((xk1(2)*delT)/2) (-(xk1(1)*delT)/2) 1];
    xin(:,i) = xk1;
    zin(:,i) = zk;
    xk = xk1;
end

%% Get estimate Kalman
%init kalman
pmin = I;
% update estimation
xestmin = xinit;
for i = 1:200
    K = pmin*H'*(H*pmin*H'+R);
    Zestmin = H*xestmin + [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001]*eye(7)*randn(7,1);
    err(:,i) = (zin(:,i)-Zestmin);
    xest = xestmin + K*err(:,i);
    zest = H*xest + [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001]*eye(7)*randn(7,1);
    P = (I-K*H)*pmin;
    Ak = [exp(-delT/tcons) zeros(1,6);
        0 exp(-delT/tcons) zeros(1,5);
        0 0 exp(-delT/tcons) zeros(1,4);
        (-(xest(5)*delT)/2) (-(xest(6)*delT)/2) (-(xest(7)*delT)/2) 1 (-(xest(1)*delT)/2) (-(xest(2)*delT)/2) (-(xest(3)*delT)/2);
        ((xest(4)*delT)/2) (-(xest(7)*delT)/2) ((xest(6)*delT)/2) ((xest(1)*delT)/2) 1 ((xest(3)*delT)/2) (-(xest(2)*delT)/2);
        ((xest(7)*delT)/2) ((xest(4)*delT)/2) (-(xest(5)*delT)/2) ((xest(2)*delT)/2) (-(xest(3)*delT)/2) 1 ((xest(1)*delT)/2);
        (-(xest(6)*delT)/2) ((xest(5)*delT)/2) ((xest(4)*delT)/2) ((xest(3)*delT)/2) ((xest(2)*delT)/2) (-(xest(1)*delT)/2) 1];
    xestmin = Ak*xest;
    pmin = Ak*P*Ak'+Q;
    Xfil(:,i) = xest;
    Zfil(:,i) = zest;
    tracePk(:,i) = trace(P);
end

%% calculate
MSE = (sum(err.^2,2)/200)
RMSE = MSE.^0.5

%% Plot Estimate vs Reference
figure(1)
title('Estimated vs Reference Roll Angular Rates Trajectory')
plot(1:200,Zfil(1,:))
hold on
plot(1:200,zin(1,:),'. -')
hold off
legend({'KALMAN','QUEST'},'FontSize',12)
xlabel('Samples (n)')
ylabel('wx (rad/sec)')
figure(2)
title('Estimated vs Reference Pitch Angular Rates Trajectory')
plot(1:200,Zfil(2,:))
hold on
plot(1:200,zin(2,:),'. -')
hold off
legend({'KALMAN','QUEST'},'FontSize',12)
xlabel('Samples (n)')
ylabel('wy (rad/sec)')
figure(3)
title('Estimated vs Reference Yaw Angular Rates Trajectory')
plot(1:200,Zfil(3,:))
hold on
plot(1:200,zin(3,:),'. -')
hold off
legend({'KALMAN','QUEST'},'FontSize',12)
xlabel('Samples (n)')
ylabel('wz (rad/sec)')
figure(4)
s(1) = subplot(2,1,1);
plot(1:200,Zfil(4:7,:))
title(s(1),'Kalman')
xlabel(s(1),'Samples (n)')
ylabel(s(1),'a,b,c,d')
s(2) = subplot(2,1,2);
plot(1:200,zin(4:7,:))
title(s(2),'QUEST')
xlabel(s(2),'Samples (n)')
ylabel(s(2),'a,b,c,d')

%% prediksi 1 (belom direkomendasikan wkwkwk)
% %init kalman
% pminpred = I;
% % update estimation
% xestminpred = xinit;
% for i = 1:200
%     Kpred = pminpred*H'*(H*pminpred*H'+R);
%     Zestminpred = H*xestminpred + [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001]*eye(7)*randn(7,1);
%     errpred(:,i) = (zin(:,i)-Zestminpred);
%     xestpred = xestminpred + Kpred*errpred(:,i);
%     zestpred = H*xestpred + [0.01 0.01 0.01 0.0001 0.0001 0.0001 0.0001]*eye(7)*randn(7,1);
%     Ppred = (I-Kpred*H)*pminpred;
%     Akpred = [exp(-delT/tcons) zeros(1,6);
%         0 exp(-delT/tcons) zeros(1,5);
%         0 0 exp(-delT/tcons) zeros(1,4);
%         (-(xestpred(5)*delT)/2) (-(xestpred(6)*delT)/2) (-(xestpred(7)*delT)/2) 1 (-(xestpred(1)*delT)/2) (-(xestpred(2)*delT)/2) (-(xestpred(3)*delT)/2);
%         ((xestpred(4)*delT)/2) (-(xestpred(7)*delT)/2) ((xestpred(6)*delT)/2) ((xestpred(1)*delT)/2) 1 ((xestpred(3)*delT)/2) (-(xestpred(2)*delT)/2);
%         ((xestpred(7)*delT)/2) ((xestpred(4)*delT)/2) (-(xestpred(5)*delT)/2) ((xestpred(2)*delT)/2) (-(xestpred(3)*delT)/2) 1 ((xestpred(1)*delT)/2);
%         (-(xestpred(6)*delT)/2) ((xestpred(5)*delT)/2) ((xestpred(4)*delT)/2) ((xestpred(3)*delT)/2) ((xestpred(2)*delT)/2) (-(xestpred(1)*delT)/2) 1];
%     xestminpred = Akpred*xestminpred;
%     pminpred = Akpred*Ppred*Akpred'+Q;
%     Xpred(:,i) = xestpred;
%     Zpred(:,i) = zestpred;
%     tracePkpred(:,i) = trace(Ppred);
% end
% 
% %% Plot Estimate vs Pred
% figure(5)
% title('Estimated vs Reference Roll Angular Rates Trajectory')
% plot(1:200,Zfil(1,:))
% hold on
% plot(1:200,Zpred(1,:),'. -')
% hold off
% legend({'KALMAN','Prediksi 1 langkah'},'FontSize',12)
% xlabel('Samples (n)')
% ylabel('wx (rad/sec)')
% figure(6)
% title('Estimated vs Reference Pitch Angular Rates Trajectory')
% plot(1:200,Zfil(2,:))
% hold on
% plot(1:200,Zpred(2,:),'. -')
% hold off
% legend({'KALMAN','Prediksi 1 langkah'},'FontSize',12)
% xlabel('Samples (n)')
% ylabel('wy (rad/sec)')
% figure(7)
% title('Estimated vs Reference Yaw Angular Rates Trajectory')
% plot(1:200,Zfil(3,:))
% hold on
% plot(1:200,Zpred(3,:),'. -')
% hold off
% legend({'KALMAN','Prediksi 1 langkah'},'FontSize',12)
% xlabel('Samples (n)')
% ylabel('wz (rad/sec)')
% figure(8)
% s(1) = subplot(2,1,1);
% plot(1:200,Zfil(4:7,:))
% title(s(1),'Kalman')
% xlabel(s(1),'Samples (n)')
% ylabel(s(1),'a,b,c,d')
% s(2) = subplot(2,1,2);
% plot(1:200,Zpred(4:7,:))
% title(s(2),'Prediksi 1 Langkah')
% xlabel(s(2),'Samples (n)')
% ylabel(s(2),'a,b,c,d')
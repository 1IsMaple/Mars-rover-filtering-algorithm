clear; close all; clc

load("CKF\error_xCKF.mat");
load("CKF\error_vCKF.mat");
load("CKF\RMSE_posCKF.mat");
load("CKF\RMSE_velCKF.mat");

load("UKF\error_xUKF.mat");
load("UKF\error_vUKF.mat");
load("UKF\RMSE_posUKF.mat");
load("UKF\RMSE_velUKF.mat");

load("EKF\error_xEKF.mat");
load("EKF\error_vEKF.mat");
load("EKF\RMSE_posEKF.mat");
load("EKF\RMSE_velEKF.mat");


%% 绘制EKF、UKF、CKF对比图
figure(1);
hold on;
grid on;
plot(error_xCKF(1:1:40),'r.-');
plot(error_xUKF(1:1:40),'g.-');
plot(error_xEKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('位置误差[km]')
title('三种滤波算法位置误差对比');
legend('CKF','UKF','EKF')
% legend('CKF','UKF')

figure(2);
hold on;
grid on;
plot(error_vCKF(1:1:40),'r.-');
plot(error_vUKF(1:1:40),'g.-');
plot(error_vEKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度误差[km/s]')
title('三种滤波算法速度误差对比');
legend('CKF','UKF','EKF')
% legend('CKF','UKF')

figure(3);
hold on;
grid on;
plot(RMSE_posCKF(1:1:40),'r.-');
plot(RMSE_posUKF(1:1:40),'g.-');
plot(RMSE_posEKF(1:1:40),'b.-'); 
xlabel('采样次数/次');
ylabel('位置RMSE[km]')
title('三种滤波算法位置RMSE对比');
legend('CKF','UKF','EKF')
% legend('CKF','UKF')

figure(4);
hold on;
grid on;
plot(RMSE_velCKF(1:1:40),'r.-');
plot(RMSE_velUKF(1:1:40),'g.-');
plot(RMSE_velEKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度RMSE[km/s]')
title('三种滤波算法速度RMSE对比');
legend('CKF','UKF','EKF')
% legend('CKF','UKF')

%%  计算仿真数据，各项数据mean
errorUKFx =0; errorUKFvx =0;RMSE_UKFx=0;RMSE_UKFvx=0;
errorCKFx=0;errorCKFvx=0;RMSE_CKFx=0;RMSE_CKFvx=0;
errorEKFx=0;errorEKFvx=0; RMSE_EKFx =0;RMSE_EKFvx =0;
for b=1:40

    errorEKFx = errorEKFx + error_xEKF(b);
    errorEKFvx= errorEKFvx + error_vEKF(b);
    RMSE_EKFx = RMSE_EKFx + RMSE_posEKF(b);
    RMSE_EKFvx = RMSE_EKFvx + RMSE_velEKF(b);

    errorUKFx = errorUKFx + error_xUKF(b);
    errorUKFvx= errorUKFvx + error_vUKF(b);
    RMSE_UKFx = RMSE_UKFx + RMSE_posUKF(b);
    RMSE_UKFvx = RMSE_UKFvx + RMSE_velUKF(b);

    errorCKFx = errorCKFx + error_xCKF(b);
    errorCKFvx= errorCKFvx + error_vCKF(b);
    RMSE_CKFx = RMSE_CKFx + RMSE_posCKF(b);
    RMSE_CKFvx = RMSE_CKFvx + RMSE_velCKF(b);
end
errorUKFx = errorUKFx/40;
errorUKFvx= errorUKFvx/40;
RMSE_UKFx = RMSE_UKFx/40;
RMSE_UKFvx = RMSE_UKFvx/40;
errorCKFx = errorCKFx/40;
errorCKFvx= errorCKFvx/40;
RMSE_CKFx = RMSE_CKFx/40;
RMSE_CKFvx = RMSE_CKFvx/40;
errorEKFx = errorEKFx/40;
errorEKFvx= errorEKFvx/40;
RMSE_EKFx = RMSE_EKFx/40;
RMSE_EKFvx = RMSE_EKFvx/40;
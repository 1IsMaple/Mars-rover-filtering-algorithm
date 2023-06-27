clear; close all; clc

load("ASRUKF\error_xASRUKF.mat");
load("ASRUKF\error_vASRUKF.mat");
load("ASRUKF\RMSE_posASRUKF.mat");
load("ASRUKF\RMSE_velASRUKF.mat");

load("UKF\error_xUKF.mat");
load("UKF\error_vUKF.mat");
load("UKF\RMSE_posUKF.mat");
load("UKF\RMSE_velUKF.mat");

%% 绘制UKF、ASRUKF对比图
figure(1);
hold on;
grid on;
plot(error_xUKF(1:1:40),'r.-');
plot(error_xASRUKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('位置误差[km]')
title('UKF、ASRUKF位置误差对比');
legend('UKF','ASRUKF')

figure(2);
hold on;
grid on;
plot(error_vUKF(1:1:40),'r.-');
plot(error_vASRUKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度误差[km/s]')
title('UKF、ASRUKF速度误差对比');
legend('UKF','ASRUKF')

figure(3);
hold on;
grid on;
plot(RMSE_posUKF(1:1:40),'r.-');
plot(RMSE_posASRUKF(1:1:40),'b.-'); 
xlabel('采样次数/次');
ylabel('位置RMSE[km]')
title('UKF、ASRUKF位置RMSE对比');
legend('UKF','ASRUKF')

figure(4);
hold on;
grid on;
plot(RMSE_velUKF(1:1:40),'r.-');
plot(RMSE_velASRUKF(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度RMSE[km/s]')
title('UKF、ASRUKF速度RMSE对比');
legend('UKF','ASRUKF')

%%  计算仿真数据，各项数据mean
errorUKFx =0; errorUKFvx =0;RMSE_UKFx=0;RMSE_UKFvx=0;
errorASRUKFx =0; errorASRUKFvx =0;RMSE_ASRUKFx=0;RMSE_ASRUKFvx=0;
for b=1:40
    errorUKFx = errorUKFx + error_xUKF(b);
    errorUKFvx= errorUKFvx + error_vUKF(b);
    RMSE_UKFx = RMSE_UKFx + RMSE_posUKF(b);
    RMSE_UKFvx = RMSE_UKFvx + RMSE_velUKF(b);

    errorASRUKFx = errorASRUKFx + error_xASRUKF(b);
    errorASRUKFvx= errorASRUKFvx + error_xASRUKF(b);
    RMSE_ASRUKFx = RMSE_ASRUKFx + RMSE_posASRUKF(b);
    RMSE_ASRUKFvx = RMSE_ASRUKFvx + RMSE_velASRUKF(b);
end
errorUKFx = errorUKFx/40;
errorUKFvx= errorUKFvx/40;
RMSE_UKFx = RMSE_UKFx/40;
RMSE_UKFvx = RMSE_UKFvx/40;
errorASRUKFx = errorASRUKFx/40;
errorASRUKFvx= errorASRUKFvx/40;
RMSE_ASRUKFx = RMSE_ASRUKFx/40;
RMSE_ASRUKFvx = RMSE_ASRUKFvx/40;

clear; close all; clc

load("ASRUKF\error_xASRUKF.mat");
load("ASRUKF\error_vASRUKF.mat");
load("ASRUKF\RMSE_posASRUKF.mat");
load("ASRUKF\RMSE_velASRUKF.mat");

load("ASRUKF001q\error_xASRUKF001.mat");
load("ASRUKF001q\error_vASRUKF001.mat");
load("ASRUKF001q\RMSE_posASRUKF001.mat");
load("ASRUKF001q\RMSE_velASRUKF001.mat");

load("ASRUKF100q\error_xASRUKF100.mat");
load("ASRUKF100q\error_vASRUKF100.mat");
load("ASRUKF100q\RMSE_posASRUKF100.mat");
load("ASRUKF100q\RMSE_velASRUKF100.mat");

%%  绘制ASRUKF在k=0.01,1,100时对比图
figure(1);
hold on;
grid on;
plot(error_xASRUKF001(1:1:40),'r.-');
plot(error_xASRUKF(1:1:40),'g.-');
plot(error_xASRUKF100(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('位置误差[km]')
title('ASRUKF在k=0.01,1,100时位置误差对比');
legend('k=0.01','k=1','k=100')

figure(2);
hold on;
grid on;
plot(error_vASRUKF001(1:1:40),'r.-');
plot(error_vASRUKF(1:1:40),'g.-');
plot(error_vASRUKF100(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度误差[km/s]')
title('ASRUKF在k=0.01,1,100时速度误差对比');
legend('k=0.01','k=1','k=100')

figure(3);
hold on;
grid on;
plot(RMSE_posASRUKF001(1:1:40),'r.-');
plot(RMSE_posASRUKF(1:1:40),'g.-'); 
plot(RMSE_posASRUKF100(1:1:40),'b.-'); 
xlabel('采样次数/次');
ylabel('位置RMSE[km]')
title('ASRUKF在k=0.01,1,100时位置RMSE');
legend('k=0.01','k=1','k=100')

figure(4);
hold on;
grid on;
plot(RMSE_velASRUKF001(1:1:40),'r.-');
plot(RMSE_velASRUKF(1:1:40),'g.-');
plot(RMSE_velASRUKF100(1:1:40),'b.-');
xlabel('采样次数/次');
ylabel('速度RMSE[km/s]')
title('ASRUKF在k=0.01,1,100时速度RMSE');
legend('k=0.01','k=1','k=100')


%%  计算仿真数据，各项数据mean
errorASRUKFx =0; errorASRUKFvx =0;RMSE_ASRUKFx=0;RMSE_ASRUKFvx=0;
errorASRUKFx001 =0; errorASRUKFvx001 =0;RMSE_ASRUKFx001=0;RMSE_ASRUKFvx001=0;
errorASRUKFx100 =0; errorASRUKFvx100 =0;RMSE_ASRUKFx100=0;RMSE_ASRUKFvx100=0;
for b=1:40

    errorASRUKFx001 = errorASRUKFx001 + error_xASRUKF001(b);
    errorASRUKFvx001= errorASRUKFvx001 + error_xASRUKF001(b);
    RMSE_ASRUKFx001 = RMSE_ASRUKFx001 + RMSE_posASRUKF001(b);
    RMSE_ASRUKFvx001 = RMSE_ASRUKFvx001 + RMSE_velASRUKF001(b);

    errorASRUKFx = errorASRUKFx + error_xASRUKF(b);
    errorASRUKFvx= errorASRUKFvx + error_xASRUKF(b);
    RMSE_ASRUKFx = RMSE_ASRUKFx + RMSE_posASRUKF(b);
    RMSE_ASRUKFvx = RMSE_ASRUKFvx + RMSE_velASRUKF(b);

    errorASRUKFx100 = errorASRUKFx100 + error_xASRUKF100(b);
    errorASRUKFvx100= errorASRUKFvx100 + error_xASRUKF100(b);
    RMSE_ASRUKFx100 = RMSE_ASRUKFx100 + RMSE_posASRUKF100(b);
    RMSE_ASRUKFvx100 = RMSE_ASRUKFvx100 + RMSE_velASRUKF100(b);
end

errorASRUKFx = errorASRUKFx/40;
errorASRUKFvx= errorASRUKFvx/40;
RMSE_ASRUKFx = RMSE_ASRUKFx/40;
RMSE_ASRUKFvx = RMSE_ASRUKFvx/40;

errorASRUKFx100 = errorASRUKFx100/40;
errorASRUKFvx100= errorASRUKFvx100/40;
RMSE_ASRUKFx100 = RMSE_ASRUKFx100/40;
RMSE_ASRUKFvx100 = RMSE_ASRUKFvx100/40;

errorASRUKFx001 = errorASRUKFx001/40;
errorASRUKFvx001= errorASRUKFvx001/40;
RMSE_ASRUKFx001 = RMSE_ASRUKFx001/40;
RMSE_ASRUKFvx001 = RMSE_ASRUKFvx001/40;
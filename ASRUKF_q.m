clear; close all; clc

filter = SingleTargetFilter;
filter = filter.gen_model;
MCRuns = 1;     %
%% 自适应开方根无迹卡尔曼滤波ASRUKF
RMSE_posASRUKFx = zeros(MCRuns,filter.K); 
RMSE_posASRUKFy = zeros(MCRuns,filter.K); 
RMSE_posASRUKFz = zeros(MCRuns,filter.K); 
RMSE_velASRUKFx = zeros(MCRuns,filter.K);
RMSE_velASRUKFy = zeros(MCRuns,filter.K);
RMSE_velASRUKFz = zeros(MCRuns,filter.K);
RMSE_posASRUKF  = zeros(MCRuns,filter.K);
RMSE_velASRUKF  = zeros(MCRuns,filter.K);

for iMCruns = 1:MCRuns
    p=[100^2;200^2;70^2;(1)^2;(1)^2;(1)^2];
    y0=[5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169];
%     y0=[10,-10,10,0.5,0.5,0.5];
%     x_pre  = [5,-5,5,0.025,0.025,0.025];
%     [5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169]
    p0=[0.002 0.0015 0.0007 0.003e-9 0.03e-9 0.03e-9];

    stateUpd_ASRUKF = (y0+p0)';
    covarUpd_ASRUKF = diag(1.*p);

    est_ASRUKF = zeros(filter.targetStateDim,filter.K);

    tASRUKF = 0;
    for k = 1:filter.K
        h = waitbar(iMCruns/MCRuns);
         %%
        tic
        % UKF预测
        [weightState_SP,statePre_ASRUFK,covarPre_ASRUKF] = filter.ASRUKFpredict(stateUpd_ASRUKF,covarUpd_ASRUKF);
        % UKF校正
        [stateUpd_ASRUKF,covarUpd_ASRUKF] = filter.ASRUKFupdate(filter.meas(:,k),statePre_ASRUFK,covarPre_ASRUKF,weightState_SP);
        % 保存滤波结果
        est_ASRUKF(:,k) = stateUpd_ASRUKF;
        RMSE_posASRUKFx(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(1)-filter.truth_X(1,k)).^2);
        RMSE_posASRUKFy(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(2)-filter.truth_X(2,k)).^2);
        RMSE_posASRUKFz(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(3)-filter.truth_X(3,k)).^2);
        RMSE_velASRUKFx(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(4)-filter.truth_X(4,k)).^2);
        RMSE_velASRUKFy(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(5)-filter.truth_X(5,k)).^2);
        RMSE_velASRUKFz(iMCruns,k) = sqrt(sum(stateUpd_ASRUKF(6)-filter.truth_X(6,k)).^2);
        RMSE_posASRUKF(iMCruns,k)=sqrt(sum(stateUpd_ASRUKF([1 2 3])-filter.truth_X([1 2 3],k)).^2);
        RMSE_velASRUKF(iMCruns,k)=sqrt(sum(stateUpd_ASRUKF([4 5 6])-filter.truth_X([4 5 6],k)).^2);

        tASRUKF = tASRUKF+toc;
    end
        disp('========================');
    disp('耗费时间/s：');
    disp(['ASRUKF:',num2str(tASRUKF)]);
    disp(k);
    disp('x位置的RMSE：');
    disp(RMSE_posASRUKFx(iMCruns,k));
    disp('y位置的RMSE：');
    disp(RMSE_posASRUKFy(iMCruns,k));
    disp('z位置的RMSE：');
    disp(RMSE_posASRUKFz(iMCruns,k));
    disp('x速度的RMSE：');
    disp(RMSE_posASRUKFx(iMCruns,k));
    disp('y速度的RMSE：');
    disp(RMSE_posASRUKFy(iMCruns,k));
    disp('z速度的RMSE：')
    disp(RMSE_posASRUKFz(iMCruns,k));
    a=sum(abs(est_ASRUKF'-filter.truth_X'));
    b=sqrt(a.^2/filter.K);
    disp(b);
end

close(h);

error_xASRU=(est_ASRUKF(1,:)-filter.truth_X(1,:));
error_yASRU=est_ASRUKF(2,:)-filter.truth_X(2,:);
error_zASRU=(est_ASRUKF(3,:)-filter.truth_X(3,:));
error_vxASRU=est_ASRUKF(4,:)-filter.truth_X(4,:);
error_vyASRU=est_ASRUKF(5,:)-filter.truth_X(5,:);
error_vzASRU=est_ASRUKF(6,:)-filter.truth_X(6,:);
sum(abs(est_ASRUKF-filter.truth_X));
figure(1)
grid on;
hold on;
plot(60:60:2400,error_xASRU(1:1:40), 'r');
plot(60:60:2400,error_yASRU(1:1:40),'g');
plot(60:60:2400,error_zASRU(1:1:40),'b');
ylabel('ASRUKF位置误差[km]')
xlabel('ASRUKF仿真时间/s')
title('ASRUKF位置误差');
legend('X','Y','Z')


figure(2)
grid on;
hold on;
plot(60:60:2400,error_vxASRU(1:1:40), 'r');
plot(60:60:2400,error_vyASRU(1:1:40),'g');
plot(60:60:2400,error_vzASRU(1:1:40),'b');

ylabel('ASRUKF速度误差[km/s]')
xlabel('ASRUKF仿真时间/s')
title('ASRUKF速度误差');
legend('X','Y','Z')
 

figure(3);
subplot(311)
plot(1:filter.K,RMSE_posASRUKFx,'r','LineWidth',1.5);
ylabel('X的RMSE[km]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的X位置RMSE');
subplot(312)
plot(1:filter.K,RMSE_posASRUKFy,'g','LineWidth',1.5);
ylabel('Y的RMSE[km]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的Y位置RMSE');
subplot(313)
plot(1:filter.K,RMSE_posASRUKFz,'b','LineWidth',1.5);
ylabel('Z的RMSE[km]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的Z位置RMSE');

figure(4);
subplot(311)
plot(1:filter.K,RMSE_velASRUKFx,'r','LineWidth',1.5);
ylabel('VX的RMSE[km/s]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的X轴速度RMSE');
subplot(312)
plot(1:filter.K,RMSE_velASRUKFy,'g','LineWidth',1.5);
ylabel('VY的RMSE[km/s]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的Y轴速度RMSE');
subplot(313)
plot(1:filter.K,RMSE_velASRUKFz,'b','LineWidth',1.5);
ylabel('VZ的RMSE[km/s]')
xlabel('ASRUKF采样次数/次')
title('ASRUKF的Z轴速度RMSE');

x=400;
for a=1:x
error_xASRUKF(a)=sqrt((error_xASRU(a))^2+(error_yASRU(a))^2+(error_zASRU(a))^2);
error_vASRUKF(a)=sqrt((error_vxASRU(a))^2+(error_vyASRU(a))^2+(error_vzASRU(a))^2);
end

save('error_xASRUKF');
save('error_vASRUKF');

save('RMSE_posASRUKF');
save('RMSE_velASRUKF');
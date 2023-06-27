clear; close all; clc

filter = SingleTargetFilter;
filter = filter.gen_model;
MCRuns = 1;     %
%% 无迹卡尔曼滤波UKF
RMSE_posUKFx = zeros(MCRuns,filter.K); 
RMSE_posUKFy = zeros(MCRuns,filter.K); 
RMSE_posUKFz = zeros(MCRuns,filter.K); 
RMSE_velUKFx = zeros(MCRuns,filter.K);
RMSE_velUKFy = zeros(MCRuns,filter.K);
RMSE_velUKFz = zeros(MCRuns,filter.K);
RMSE_posUKF  = zeros(MCRuns,filter.K);
RMSE_velUKF  = zeros(MCRuns,filter.K);

for iMCruns = 1:MCRuns
    p=[100^2;200^2;70^2;(1)^2;(1)^2;(1)^2];
    y0=[5,-5,5,0.025,0.025,0.025];
%     y0=[10,-10,10,0.5,0.5,0.5];
%     x_pre  = [5,-5,5,0.025,0.025,0.025];
%     [5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169]
    p0=[0.002 0.0015 0.0007 0.003e-9 0.03e-9 0.03e-9];

    stateUpd_UKF = (y0+p0)';
    covarUpd_UKF = diag(1.*p);

    est_UKF = zeros(filter.targetStateDim,filter.K);

    tUKF = 0;
    for k = 1:filter.K
        h = waitbar(iMCruns/MCRuns);
         %%
        tic
        % UKF预测
        [weightState_SP,statePre_UFK,covarPre_UKF] = filter.UKFpredict(stateUpd_UKF,covarUpd_UKF);
        % UKF校正
        [stateUpd_UKF,covarUpd_UKF] = filter.UKFupdate(filter.meas(:,k),statePre_UFK,covarPre_UKF,weightState_SP);
        % 保存滤波结果
        est_UKF(:,k) = stateUpd_UKF;
        RMSE_posUKFx(iMCruns,k) = sqrt(sum(stateUpd_UKF(1)-filter.truth_X(1,k)).^2);
        RMSE_posUKFy(iMCruns,k) = sqrt(sum(stateUpd_UKF(2)-filter.truth_X(2,k)).^2);
        RMSE_posUKFz(iMCruns,k) = sqrt(sum(stateUpd_UKF(3)-filter.truth_X(3,k)).^2);
        RMSE_velUKFx(iMCruns,k) = sqrt(sum(stateUpd_UKF(4)-filter.truth_X(4,k)).^2);
        RMSE_velUKFy(iMCruns,k) = sqrt(sum(stateUpd_UKF(5)-filter.truth_X(5,k)).^2);
        RMSE_velUKFz(iMCruns,k) = sqrt(sum(stateUpd_UKF(6)-filter.truth_X(6,k)).^2);
        RMSE_posUKF(iMCruns,k)=sqrt(sum(stateUpd_UKF([1 2 3])-filter.truth_X([1 2 3],k)).^2);
        RMSE_velUKF(iMCruns,k)=sqrt(sum(stateUpd_UKF([4 5 6])-filter.truth_X([4 5 6],k)).^2);
        tUKF = tUKF+toc;
    end
        disp('========================');
    disp('耗费时间/s：');
    disp(['UKF:',num2str(tUKF)]);
    disp(k);
    disp('x位置的RMSE：');
    disp(RMSE_posUKFx(iMCruns,k));
    disp('y位置的RMSE：');
    disp(RMSE_posUKFy(iMCruns,k));
    disp('z位置的RMSE：');
    disp(RMSE_posUKFz(iMCruns,k));
    disp('x速度的RMSE：');
    disp(RMSE_posUKFx(iMCruns,k));
    disp('y速度的RMSE：');
    disp(RMSE_posUKFy(iMCruns,k));
    disp('z速度的RMSE：')
    disp(RMSE_posUKFz(iMCruns,k));
    a=sum(abs(est_UKF'-filter.truth_X'));
    b=sqrt(a.^2/filter.K);
    disp(b);
end

close(h);

error_xU=(est_UKF(1,:)-filter.truth_X(1,:));
error_yU=est_UKF(2,:)-filter.truth_X(2,:);
error_zU=(est_UKF(3,:)-filter.truth_X(3,:));
error_vxU=est_UKF(4,:)-filter.truth_X(4,:);
error_vyU=est_UKF(5,:)-filter.truth_X(5,:);
error_vzU=est_UKF(6,:)-filter.truth_X(6,:);
sum(abs(est_UKF-filter.truth_X));
figure(1)
grid on;
hold on;
plot(60:60:24000,error_xU(1:1:end), 'r');
plot(60:60:24000,error_yU(1:1:end),'g');
plot(60:60:24000,error_zU(1:1:end),'b');
ylabel('UKF位置误差[km]')
xlabel('UKF仿真时间/s')
title('UKF位置误差');
legend('X','Y','Z')
axes('Position',[0.38,0.6,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_xU(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_yU(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_zU(1:1:end),'b');
xlim([15000,22000]); % 设置坐标轴范围  
figure(2)
grid on;
hold on;
plot(60:60:24000,error_vxU(1:1:end), 'r');
plot(60:60:24000,error_vyU(1:1:end),'g');
plot(60:60:24000,error_vzU(1:1:end),'b');
ylabel('UKF速度误差[km/s]')
xlabel('UKF仿真时间/s')
title('UKF速度误差');
legend('X','Y','Z')
axes('Position',[0.4,0.5,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_vxU(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_vyU(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_vzU(1:1:end),'b');   
xlim([15000,22000]); % 设置坐标轴范围  

figure(3);
subplot(311)
plot(1:filter.K,RMSE_posUKFx,'r','LineWidth',1.5);
ylabel('X的RMSE[km]')
xlabel('UKF采样次数/次')
title('UKF的X位置RMSE');
subplot(312)
plot(1:filter.K,RMSE_posUKFy,'g','LineWidth',1.5);
ylabel('Y的RMSE[km]')
xlabel('UKF采样次数/次')
title('UKF的Y位置RMSE');
subplot(313)
plot(1:filter.K,RMSE_posUKFz,'b','LineWidth',1.5);
ylabel('Z的RMSE[km]')
xlabel('UKF采样次数/次')
title('UKF的Z位置RMSE');

figure(4);
subplot(311)
plot(1:filter.K,RMSE_velUKFx,'r','LineWidth',1.5);
ylabel('VX的RMSE[km/s]')
xlabel('UKF采样次数/次')
title('UKF的X轴速度RMSE');
subplot(312)
plot(1:filter.K,RMSE_velUKFy,'g','LineWidth',1.5);
ylabel('VY的RMSE[km/s]')
xlabel('UKF采样次数/次')
title('UKF的Y轴速度RMSE');
subplot(313)
plot(1:filter.K,RMSE_velUKFz,'b','LineWidth',1.5);
ylabel('VZ的RMSE[km/s]')
xlabel('UKF采样次数/次')
title('UKF的Z轴速度RMSE');

x=400;
for a=1:x
error_xUKF(a)=sqrt((error_xU(a))^2+(error_yU(a))^2+(error_zU(a))^2);
error_vUKF(a)=sqrt((error_vxU(a))^2+(error_vyU(a))^2+(error_vzU(a))^2);
end

save('error_xUKF');
save('error_vUKF');

save('RMSE_posUKF');
save('RMSE_velUKF');
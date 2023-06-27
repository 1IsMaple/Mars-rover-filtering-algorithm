clear; close all; clc

filter = SingleTargetFilter;
filter = filter.gen_model;
MCRuns = 1;     %
%% 扩展卡尔曼滤波EKF
RMSE_posEKFx = zeros(MCRuns,filter.K);
RMSE_posEKFy = zeros(MCRuns,filter.K);
RMSE_posEKFz = zeros(MCRuns,filter.K);
RMSE_velEKFx = zeros(MCRuns,filter.K);
RMSE_velEKFy = zeros(MCRuns,filter.K);
RMSE_velEKFz = zeros(MCRuns,filter.K);
RMSE_posEKF  = zeros(MCRuns,filter.K);
RMSE_velEKF  = zeros(MCRuns,filter.K);


for iMCruns = 1:MCRuns
    p=[100^2;200^2;70^2;(1)^2;(1)^2;(1)^2];
    y0=[5,-5,5,0.025,0.025,0.025];
%   y0=[10,-10,10,0.5,0.5,0.5];
%   x_pre  = [5,-5,5,0.025,0.025,0.025];
%   [5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169]
    p0=[0.002 0.0015 0.0007 0.003e-9 0.03e-9 0.03e-9];

    stateUpd_EKF = (y0+p0)';
    covarUpd_EKF = diag(1.*p);

    est_EKF = zeros(filter.targetStateDim,filter.K);

    tEKF = 0; tUKF = 0; tCKF = 0; tPF = 0;
    for k = 1:filter.K
        h = waitbar(iMCruns/MCRuns);
        %%
        tic
        % EKF预测
        [statePre_EFK,covarPre_EKF] = filter.EKFpredict(stateUpd_EKF,covarUpd_EKF);
        % EKF校正
        [stateUpd_EKF,covarUpd_EKF] = filter.EKFupdate(filter.meas(:,k),statePre_EFK,covarPre_EKF);
        % 保存滤波结果
        est_EKF(:,k) = stateUpd_EKF;
        RMSE_posEKFx(iMCruns,k) = sqrt(sum(stateUpd_EKF(1)-filter.truth_X(1,k)).^2);
        RMSE_posEKFy(iMCruns,k) = sqrt(sum(stateUpd_EKF(2)-filter.truth_X(2,k)).^2);
        RMSE_posEKFz(iMCruns,k) = sqrt(sum(stateUpd_EKF(3)-filter.truth_X(3,k)).^2);

        RMSE_velEKFx(iMCruns,k) = sqrt(sum(stateUpd_EKF(4)-filter.truth_X(4,k)).^2);
        RMSE_velEKFy(iMCruns,k) = sqrt(sum(stateUpd_EKF(5)-filter.truth_X(5,k)).^2);
        RMSE_velEKFz(iMCruns,k) = sqrt(sum(stateUpd_EKF(6)-filter.truth_X(6,k)).^2);

        RMSE_posEKF(iMCruns,k)=sqrt(sum(stateUpd_EKF([1 2 3])-filter.truth_X([1 2 3],k)).^2);
        RMSE_velEKF(iMCruns,k)=sqrt(sum(stateUpd_EKF([4 5 6])-filter.truth_X([4 5 6],k)).^2);

        tEKF = tEKF+toc;

    end
        disp('========================');
    disp('耗费时间/s：');
    disp(['EKF:',num2str(tEKF)]);

    disp(k);
    disp('x位置的RMSE：');
    disp(RMSE_posEKFx(iMCruns,k));
    disp('y位置的RMSE：');
    disp(RMSE_posEKFy(iMCruns,k));
    disp('z位置的RMSE：');
    disp(RMSE_posEKFz(iMCruns,k));
    disp('x速度的RMSE：');
    disp(RMSE_posEKFx(iMCruns,k));
    disp('y速度的RMSE：');
    disp(RMSE_posEKFy(iMCruns,k));
    disp('z速度的RMSE：')
    disp(RMSE_posEKFz(iMCruns,k));
    a=sum(abs(est_EKF'-filter.truth_X'));
%     b=sqrt(a.^2/filter.K);   
    disp(a);
end

close(h);

error_xE=(est_EKF(1,:)-filter.truth_X(1,:));
error_yE=est_EKF(2,:)-filter.truth_X(2,:);
error_zE=(est_EKF(3,:)-filter.truth_X(3,:));
error_vxE=est_EKF(4,:)-filter.truth_X(4,:);
error_vyE=est_EKF(5,:)-filter.truth_X(5,:);
error_vzE=est_EKF(6,:)-filter.truth_X(6,:);
sum(abs(est_EKF-filter.truth_X));
figure(1)
grid on;
hold on;
plot(60:60:24000,error_xE(1:1:end), 'r');
plot(60:60:24000,error_yE(1:1:end),'g');
plot(60:60:24000,error_zE(1:1:end),'b');
ylabel('EKF位置误差[km]')
xlabel('EKF仿真时间/s')
title('EKF位置误差');
legend('X','Y','Z')
axes('Position',[0.38,0.6,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_xE(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_yE(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_zE(1:1:end),'b');
xlim([15000,22000]); % 设置坐标轴范围  
figure(2)
grid on;
hold on;
plot(60:60:24000,error_vxE(1:1:end), 'r');
plot(60:60:24000,error_vyE(1:1:end),'g');
plot(60:60:24000,error_vzE(1:1:end),'b');
ylabel('EKF速度误差[km/s]')
xlabel('EKF仿真时间/s')
title('EKF速度误差');
legend('X','Y','Z')
axes('Position',[0.4,0.5,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_vxE(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_vyE(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_vzE(1:1:end),'b');   
xlim([15000,22000]); % 设置坐标轴范围  


figure(3);
subplot(311)
plot(1:filter.K,RMSE_posEKFx,'r','LineWidth',1.5);
ylabel('X的RMSE[km]')
xlabel('EKF采样次数/次')
title('EKF的X位置RMSE');
subplot(312)
plot(1:filter.K,RMSE_posEKFy,'g','LineWidth',1.5);
ylabel('Y的RMSE[km]')
xlabel('EKF采样次数/次')
title('EKF的Y位置RMSE');
subplot(313)
plot(1:filter.K,RMSE_posEKFz,'b','LineWidth',1.5);
ylabel('Z的RMSE[km]')
xlabel('EKF采样次数/次')
title('EKF的Z位置RMSE');

figure(4);
subplot(311)
plot(1:filter.K,RMSE_velEKFx,'r','LineWidth',1.5);
ylabel('VX的RMSE[km/s]')
xlabel('EKF采样次数/次')
title('EKF的X轴速度RMSE');
subplot(312)
plot(1:filter.K,RMSE_velEKFy,'g','LineWidth',1.5);
ylabel('VY的RMSE[km/s]')
xlabel('EKF采样次数/次')
title('EKF的Y轴速度RMSE');
subplot(313)
plot(1:filter.K,RMSE_velEKFz,'b','LineWidth',1.5);
ylabel('VZ的RMSE[km/s]')
xlabel('EKF采样次数/次')
title('EKF的Z轴速度RMSE');

x=400;
for a=1:x
error_xEKF(a)=sqrt((error_xE(a))^2+(error_yE(a))^2+(error_zE(a))^2);
error_vEKF(a)=sqrt((error_vxE(a))^2+(error_vyE(a))^2+(error_vzE(a))^2);
end

save('error_xEKF');
save('error_vEKF');

save('RMSE_posEKF');
save('RMSE_velEKF');
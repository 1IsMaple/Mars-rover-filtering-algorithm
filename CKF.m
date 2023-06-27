clear; close all; clc

filter = SingleTargetFilter;
filter = filter.gen_model;
MCRuns = 1;     %
%% 容积卡尔曼滤波CKF

RMSE_posCKFx = zeros(MCRuns,filter.K); 
RMSE_posCKFy = zeros(MCRuns,filter.K); 
RMSE_posCKFz = zeros(MCRuns,filter.K); 
RMSE_velCKFx = zeros(MCRuns,filter.K);
RMSE_velCKFy = zeros(MCRuns,filter.K);
RMSE_velCKFz = zeros(MCRuns,filter.K);
RMSE_posCKF  = zeros(MCRuns,filter.K);
RMSE_velCKF  = zeros(MCRuns,filter.K);

for iMCruns = 1:MCRuns
    p=[100^2;200^2;70^2;(1)^2;(1)^2;(1)^2];
    y0=[5,-5,5,0.025,0.025,0.025];
%     y0=[10,-10,10,0.5,0.5,0.5];
    %          x_pre  = [5,-5,5,0.025,0.025,0.025];
%     [5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169]
    p0=[0.002 0.0015 0.0007 0.003e-9 0.03e-9 0.03e-9];

    stateUpd_CKF = (y0+p0)';
    covarUpd_CKF = diag(1.*p);

    est_CKF = zeros(filter.targetStateDim,filter.K);

    tCKF = 0; 
    for k = 1:filter.K
        h = waitbar(iMCruns/MCRuns);
                %%
        tic
        % CKF预测
        [statePre_CFK,covarPre_CKF] = filter.CKFpredict(stateUpd_CKF,covarUpd_CKF);
        % CKF校正
        [stateUpd_CKF,covarUpd_CKF] = filter.CKFupdate(filter.meas(:,k),statePre_CFK,covarPre_CKF);
        % 保存滤波结果
        est_CKF(:,k) = stateUpd_CKF;
        RMSE_posCKFx(iMCruns,k) = sqrt(sum(stateUpd_CKF(1)-filter.truth_X(1,k)).^2);
        RMSE_posCKFy(iMCruns,k) = sqrt(sum(stateUpd_CKF(2)-filter.truth_X(2,k)).^2);
        RMSE_posCKFz(iMCruns,k) = sqrt(sum(stateUpd_CKF(3)-filter.truth_X(3,k)).^2);
        RMSE_velCKFx(iMCruns,k) = sqrt(sum(stateUpd_CKF(4)-filter.truth_X(4,k)).^2);
        RMSE_velCKFy(iMCruns,k) = sqrt(sum(stateUpd_CKF(5)-filter.truth_X(5,k)).^2);
        RMSE_velCKFz(iMCruns,k) = sqrt(sum(stateUpd_CKF(6)-filter.truth_X(6,k)).^2);
        RMSE_posCKF(iMCruns,k)=sqrt(sum(stateUpd_CKF([1 2 3])-filter.truth_X([1 2 3],k)).^2);
        RMSE_velCKF(iMCruns,k)=sqrt(sum(stateUpd_CKF([4 5 6])-filter.truth_X([4 5 6],k)).^2);
        tCKF = tCKF+toc;
    end
        disp('========================');
    disp('耗费时间/s：');
    disp(['CKF:',num2str(tCKF)]);
    disp('x位置的RMSE：');
    disp(RMSE_posCKFx(iMCruns,k));
    disp('y位置的RMSE：');
    disp(RMSE_posCKFy(iMCruns,k));
    disp('z位置的RMSE：');
    disp(RMSE_posCKFz(iMCruns,k));
    disp('x速度的RMSE：');
    disp(RMSE_posCKFx(iMCruns,k));
    disp('y速度的RMSE：');
    disp(RMSE_posCKFy(iMCruns,k));
    disp('z速度的RMSE：')
    disp(RMSE_posCKFz(iMCruns,k));
    a=sum(abs(est_CKF'-filter.truth_X'));
    b=sqrt(a.^2/filter.K);
    disp(b);
end

close(h);


error_xC=(est_CKF(1,:)-filter.truth_X(1,:));
error_yC=est_CKF(2,:)-filter.truth_X(2,:);
error_zC=(est_CKF(3,:)-filter.truth_X(3,:));
error_vxC=est_CKF(4,:)-filter.truth_X(4,:);
error_vyC=est_CKF(5,:)-filter.truth_X(5,:);
error_vzC=est_CKF(6,:)-filter.truth_X(6,:);

figure(1)
grid on;
hold on;
plot(60:60:24000,error_xC(1:1:end), 'r');
plot(60:60:24000,error_yC(1:1:end),'g');
plot(60:60:24000,error_zC(1:1:end),'b');
ylabel('CKF位置误差[km]')
xlabel('CKF仿真时间/s')
title('CKF位置误差');
legend('X','Y','Z')
axes('Position',[0.38,0.6,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_xC(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_yC(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_zC(1:1:end),'b');
xlim([15000,22000]); % 设置坐标轴范围  
figure(2)
grid on;
hold on;
plot(60:60:24000,error_vxC(1:1:end), 'r');
plot(60:60:24000,error_vyC(1:1:end),'g');
plot(60:60:24000,error_vzC(1:1:end),'b');
ylabel('CKF速度误差[km/s]')
xlabel('CKF仿真时间/s')
title('CKF速度误差');
legend('X','Y','Z')
axes('Position',[0.4,0.5,0.38,0.28]); % 生成子图   左右  上下 宽窄
hold on;
plot(60:60:24000,error_vxC(1:1:end),'r');                                                                                                         
plot(60:60:24000,error_vyC(1:1:end),'g');                                                                                                      
plot(60:60:24000,error_vzC(1:1:end),'b');   
xlim([15000,22000]); % 设置坐标轴范围  

% figure;
% plot(1:filter.K,RMSE_posCKF,'b.-','LineWidth',1.5);
% title('CKF位置误差');
% figure;
% plot(1:filter.K,RMSE_velCKF,'b.-','LineWidth',1.5);
% title('CKF速度误差');

figure(3);
subplot(311)
plot(1:filter.K,RMSE_posCKFx,'r','LineWidth',1.5);
ylabel('X的RMSE[km]')
xlabel('CKF采样次数/次')
title('CKF的X位置RMSE');
subplot(312)
plot(1:filter.K,RMSE_posCKFy,'g','LineWidth',1.5);
ylabel('Y的RMSE[km]')
xlabel('CKF采样次数/次')
title('CKF的Y位置RMSE');
subplot(313)
plot(1:filter.K,RMSE_posCKFz,'b','LineWidth',1.5);
ylabel('Z的RMSE[km]')
xlabel('CKF采样次数/次')
title('CKF的Z位置RMSE');

figure(4);
subplot(311)
plot(1:filter.K,RMSE_velCKFx,'r','LineWidth',1.5);
ylabel('VX的RMSE[km/s]')
xlabel('CKF采样次数/次')
title('CKF的X轴速度RMSE');
subplot(312)
plot(1:filter.K,RMSE_velCKFy,'g','LineWidth',1.5);
ylabel('VY的RMSE[km/s]')
xlabel('CKF采样次数/次')
title('CKF的Y轴速度RMSE');
subplot(313)
plot(1:filter.K,RMSE_velCKFz,'b','LineWidth',1.5);
ylabel('VZ的RMSE[km/s]')
xlabel('CKF采样次数/次')
title('CKF的Z轴速度RMSE');
x=400;
for a=1:x
error_xCKF(a)=sqrt((error_xC(a))^2+(error_yC(a))^2+(error_zC(a))^2);
error_vCKF(a)=sqrt((error_vxC(a))^2+(error_vyC(a))^2+(error_vzC(a))^2);
end

save('error_xCKF');
save('error_vCKF');

save('RMSE_posCKF');
save('RMSE_velCKF');


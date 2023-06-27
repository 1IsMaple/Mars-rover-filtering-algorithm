classdef SingleTargetFilter
    properties
        K; % 观测总帧数
        T; % 滤波间隔
        % ====== 动态模型参数 =======
        targetStateDim; % 状态维数
        sigma_process;  % 过程噪声协方差   
        processNoiseDim;% 过程噪声维数
        Q;              % 过程噪声驱动矩阵
        % ====== 测量模型参数 =======
        R;              % 测量噪声协方差
        MeasNoiseDim;   % 测量噪声维数
        MeasDim;        % 测量数据维数
        % ====== 交互参数 =======
        truth_X;        % 真实航迹
        meas;           % 观测数据
        % ====== UKF参数 ========
        alpha = 0.1;   % 0 <  alpha <= 1     (1e-3)
        beta = 2;       % 0 <= beta           (2)
        kappa = 0;      % 0 <= kappa <= 3     (0)
    end
    %% 方法
    methods
        function obj = SingleTargetFilter
            % 构造函数
            % ============== 参数初始化 ===============
            q1=1e-2;
            q21=1e-2;q22=1e-2;q23=1e-2;
            q=[q1^2;q1^2;q1^2;q21^2;q22^2;q23^2];
%           rr1=[ones(1,3).*2e-2 3e-3 3e-3 3e-3 3e-3 3e-3];
%           rr1=[ones(1,3).*2e-1 3e-1 3e-1 3e-1 3e-1 3e-1];
            rr1=[ones(1,3).*2e-2 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3];
            obj.K = 400;
            obj.T = 1;
            obj.targetStateDim = 6;
            obj.processNoiseDim = 6;
            obj.MeasNoiseDim = 11;
            obj.MeasDim = 11;
            obj.sigma_process = diag(q.^2);
            obj.R = diag(rr1.^2).*1;
            obj.Q = diag(q.^2);
        end
        function obj = gen_model(obj)          
            x_pre=[5,-5,5,0.025,0.025,0.025]';
%           x_pre  = [5,-5,5,0.025,0.025,0.025];
%           [5767802.3,-2818455,5876662.9,-4119.834,-572.527,-1294.169]
%           x_pre=[10,-10,10,0.5,0.5,0.5]';
            obj.truth_X(:,1)=x_pre;
            obj.meas(:,1)=[0;0;0;0;0;0;0;0;0;0;0];
            for k =1:obj.K
                [X_k,Z_k,Truth] = obj.f(x_pre);
                % 测量
                obj.meas(:,k) = Z_k;
                %预测：F*xk-1
                x_pre = X_k(x_pre);
                obj.truth_X(:,k) = Truth;
            end
             
        end
                  %% ============== 无迹卡尔曼滤波UKF ===============
        % 生成sigma点
        function [weightState_SP,state_SP] = gen_sigmaPoint(obj,statePrior_SP,covarPrior_SP)     
            lambda = obj.alpha^2*(obj.targetStateDim+obj.kappa)-obj.targetStateDim;
            % 协方差矩阵方根
            cholesky_decom = obj.cholPSD(covarPrior_SP*(obj.targetStateDim+lambda));
            % sigma点状态
            state_SP = zeros(obj.targetStateDim,2*obj.targetStateDim);
            for iSigmaPoint = 1:obj.targetStateDim
                state_SP(:,iSigmaPoint) = statePrior_SP+cholesky_decom(:,iSigmaPoint);
                state_SP(:,iSigmaPoint+obj.targetStateDim) = statePrior_SP-cholesky_decom(:,iSigmaPoint);
            end
            state_SP = cat(2,statePrior_SP,state_SP);
            % sigma点权重
            weightState_SP(1,:) = cat(2,lambda/(obj.targetStateDim+lambda),1/(2*(obj.targetStateDim+lambda)).*ones(1,2*obj.targetStateDim));
            weightState_SP(2,:) = cat(2,lambda/(obj.targetStateDim+lambda)+(1-obj.alpha^2+obj.beta),1/(2*(obj.targetStateDim+lambda)).*ones(1,2*obj.targetStateDim));
        end
        % 预测
        function [weight_SP,statePre,covarPre] = UKFpredict(obj,statePrior,covarPrior)
            % 获取先验状态的sigma点
            [weight_SP,state_SP] = obj.gen_sigmaPoint(statePrior,covarPrior);
            % sigma点状态预测
            sigmaPoint_num = size(state_SP,2);
            statePre_SP = zeros(obj.targetStateDim,sigmaPoint_num);
            % 逐点预测
            for i = 1:sigmaPoint_num
                statePrior = state_SP(:,i);
                statePre_SP(:,i) = obj.Mars_model(statePrior)*statePrior;
            end
            % 加权求和
            statePre = statePre_SP*weight_SP(1,:)';
            % 预测协方差
            covarPre = zeros(obj.targetStateDim,obj.targetStateDim);
            for i = 1:sigmaPoint_num
                covarPre = covarPre+weight_SP(2,i)*(statePre-statePre_SP(:,i))*(statePre-statePre_SP(:,i))';
            end
            covarPre = covarPre+obj.Q;
%             state_weightedBias = (statePre_SP-statePre).*(ones(obj.targetStateDim,1)*sqrt(weightCovar_SP));
%             covarPre = state_weightedBias*state_weightedBias'+obj.Q*sqrtm(obj.sigma_process)*obj.Q';
        end
        % 更新
        function [stateUpd,covarUpd] = UKFupdate(obj,measZ,statePre,covarPre,weight_SP)
            % 获取先验状态的sigma点
            [~,state_SP] = obj.gen_sigmaPoint(statePre,covarPre);
            % 观测预测
%             measPre_SP = [atan2(state_SP(3,:),state_SP(1,:)); sqrt(sum(state_SP([1 3],:).^2,1))];
            Ra = 2.27e8;    % 火星轨道半径
            c=1e7;
            h=@(y)[y(1) y(2) y(3) y(4) y(5) y(6) y(1)/y(3) y(2)/y(3) (4*y(1)^2+4*y(2)^2+4*y(3)^2)^(0.5) (y(1)^2+y(2)^2+y(3)^2)^(0.5) ((y(1)+Ra)^2+y(2)^2+y(3)^2)^(0.5)/c];
            sigmaPoint_num = size(state_SP,2); 
            measPre_SP=zeros(obj.MeasDim,sigmaPoint_num);
            for i=1:sigmaPoint_num        
            measPre_SP(:,i) = [h(state_SP(:,i))];
            end
            % 观测均值
            measPre = measPre_SP*weight_SP(1,:)';
            % 观测预测协方差
            covarMeas = zeros(obj.MeasDim,obj.MeasDim);
            sigmaPoint_num = size(state_SP,2);
            for i = 1:sigmaPoint_num
                covarMeas = covarMeas+weight_SP(2,i)*(measPre-measPre_SP(:,i))*(measPre-measPre_SP(:,i))';
            end
            covarMeas = covarMeas+obj.R;

            covar_StateMeas = zeros(obj.targetStateDim,obj.MeasDim);
            for i = 1:sigmaPoint_num
                covar_StateMeas = covar_StateMeas+weight_SP(2,i)*(statePre-state_SP(:,i))*(measPre-measPre_SP(:,i))';
            end
            % 卡尔曼增益
            K_gain = covar_StateMeas*pinv(covarMeas+eye(obj.MeasDim)*(1e-3)); %#ok<MINV> 
            % 状态更新
            stateUpd = statePre+K_gain*(measZ-measPre);
            % 协方差更新
            covarUpd = covarPre-K_gain*covarMeas*K_gain';
        end

    end
             %% 动态模型
    methods(Static)
        %% ============== 火星探测器近阶段模型（Mars_model） ===================
        function F = Mars_model(x_pre)
            %====Input====
            %x_pre:先验状态
            %F：根据状态方程获得的Jacobian矩阵
            ra=2.27e8;%火星轨道半径
            ra_dot =9e-9;%径向速度
            f_dot =1.05859e-7;%火星真近点角速度
            u=1.327e11;%太阳引力系数 1.327*10^11 km^3/s^2
            x=x_pre(1);
            y=x_pre(2);
            z=x_pre(3);
            F =[0,0,0,1,0,0;
                0,0,0,0,1,0;
                0,0,0,0,0,1;
                128550*x^2/(x^2+y^2+z^2)^(5/2)-42850/(x^2+y^2+z^2)^(3/2)+f_dot^2 - u/((ra + x)^2 + y^2 + z^2)^(3/2) + (3*u*(2*ra + 2*x)*(ra + x))/(2*((ra + x)^2 + y^2 + z^2)^(5/2)),128550*x*y/(x^2+y^2+z^2)^(5/2)+(3*u*y*(ra + x))/((ra + x)^2 + y^2 + z^2)^(5/2) - (2*f_dot*ra_dot)/ra,128550*x*z/(x^2+y^2+z^2)/(x^2+y^2+z^2)^(5/2)+(3*u*z*(ra + x))/((ra + x)^2 + y^2 + z^2)^(5/2),0,2*f_dot,0;
                128550*x*y/(x^2+y^2+z^2)^(5/2)+(2*f_dot*ra_dot)/ra + (3*u*y*(2*ra + 2*x))/(2*((ra + x)^2 + y^2 + z^2)^(5/2)),128550*y^2/(x^2+y^2+z^2)^(5/2)-42850/(x^2+y^2+z^2)^(3/2)+f_dot^2 - u/((ra + x)^2 + y^2 + z^2)^(3/2) + (3*u*y^2)/((ra + x)^2 + y^2 + z^2)^(5/2),128550*z*y/(x^2+y^2+z^2)^(5/2)+(3*u*y*z)/((ra + x)^2 + y^2 + z^2)^(5/2),-2*f_dot,0,0;
                128550*x*z/(x^2+y^2+z^2)^(5/2)+(3*u*z*(2*ra + 2*x))/(2*((ra + x)^2 + y^2 + z^2)^(5/2)),128550*y*z/(x^2+y^2+z^2)^(5/2)+(3*u*y*z)/((ra + x)^2 + y^2 + z^2)^(5/2),128550*z^2/(x^2+y^2+z^2)^(5/2)-42850/(x^2+y^2+z^2)^(3/2)-(u*(ra^2 + 2*ra*x + x^2 + y^2 - 2*z^2))/((ra + x)^2 + y^2 + z^2)^(5/2),0,0,0];     
            F=F+eye(6,6);
        end
        function H = Mars_measure(x_pre)
            %===Input===
            %x_pre:先验状态
            %H：观测方程在𝑘 − 1时刻的 Jacobian 矩阵。
            ra=2.27e8;%火星轨道半径
            c=1e7;
            x=x_pre(1);
            y=x_pre(2);
            z=x_pre(3);
            H =zeros(11,6);
            H(1,1) = 1;
            H(2,2) = 1;
            H(3,3) = 1;
            H(4,4) = 1;
            H(5,5) = 1;
            H(6,6) = 1;
            H(7,1) = 1/z;
            H(7,3) = -x/z^2;
            H(8,2) = 1/z;
            H(8,3) = -y/z^2;
            H(9,1) = 4*x/(4*x^2+4*y^2+4*z^2)^(0.5);
            H(9,2) = 4*y/(4*x^2+4*y^2+4*z^2)^(0.5);
            H(9,3) = 4*z/(4*x^2+4*y^2+4*z^2)^(0.5);
            H(10,1) = x/(x^2+y^2+z^2)^(0.5);
            H(10,2) = y/(x^2+y^2+z^2)^(0.5);
            H(10,3) = z/(x^2+y^2+z^2)^(0.5);
            H(11,1) = (x+ra)/(3e5*((x+ra)^2+y^2+z^2)^(0.5));
            H(11,2) = y/(3e5*((x+ra)^2+y^2+z^2)^(0.5));
            H(11,3) = z/(3e5*((x+ra)^2+y^2+z^2)^(0.5));
        end
            function [X_k,Z_k,Truth] = f(x_pre)
           % %===求Xk和Zk===
            %参数
            a = 1.05859e-7; % 小行星真近点角速率
            Va = 9e-9;      % 径向速度
            Ra = 2.27e8;    % 火星轨道半径
            u = 1.327e11;   % 太阳的万有引力常数
            u1 = 4.285e4;   % 火星的万有引力常数
            c=1e7;
            h=1;
            dy=@(y)[y(4),y(5),y(6),-u1*y(1)/(sqrt(y(1)^2+y(2)^2+y(3)^2)^3)+2*a*(y(5)-y(2)*Va/Ra)+y(1)*a^2+u/Ra^2-u*(Ra+y(1))/(sqrt((Ra+y(1))^2+y(2)^2+y(3)^2))^3,-u1*y(2)/(sqrt(y(1)^2+y(2)^2+y(3)^2)^3)-2*a*(y(4)-y(1)*Va/Ra)+y(2)*a^2-u*y(2)/(sqrt((Ra+y(1))^2+y(2)^2+y(3)^2))^3,-u1*y(3)/(sqrt(y(1)^2+y(2)^2+y(3)^2)^3)-u*y(3)/(sqrt((Ra+y(1))^2+y(2)^2+y(3)^2))^3];
            X_k=@(y)y+h*dy(y)';
            h1=@(y)[y(1) y(2) y(3) y(4) y(5) y(6) y(1)/y(3) y(2)/y(3) (4*y(1)^2+4*y(2)^2+4*y(3)^2)^(0.5) (y(1)^2+y(2)^2+y(3)^2)^(0.5) ((y(1)+Ra)^2+y(2)^2+y(3)^2)^(0.5)/c]';
            k1=dy(x_pre)';
            q1=1e-2;
            q21=1e-2;q22=1e-2;q23=1e-2;
            q=[q1^2;q1^2;q1^2;q21^2;q22^2;q23^2];
            rr1=[ones(1,3).*2e-2 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3 3e-3];
            Truth=x_pre+k1+q.*randn(6,1);
            Z_k=h1(Truth)+(rr1)'.*randn(11,1);
            end
                                % ============== 矩阵cholesky分解 ==============
        function Xi = cholPSD(A)
            [~, flag] = chol(A);
            if (flag == 0)
                Xi = (chol(A)).';
            else
                [~,S,V] = svd(A);
                Ss = sqrt(S);
                Xi = V*Ss;
            end
        end
    end
end
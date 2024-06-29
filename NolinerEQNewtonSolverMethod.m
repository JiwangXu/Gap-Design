function [ex1 , ey1 , ex2 , ey2 , Fsw , iter , error , h , pOil , tautheta , tauz , pContact , PFriction] =...
    NolinerEQNewtonSolverMethod( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
    pHigh , pLow , phi_rad , tolerance , uInitial)
%求任意位置的柱塞姿态和斜盘力
%使用切线法（牛顿法）迭代求解，相比割线法，需要使用偏导数计算雅克比矩阵，时间复杂度高，但迭代步数会减少；同时切线法更容易寻找初始解；
%输入参数：
%phi_rad：当前工作角度
%h0：平均油膜厚度；
%Fhtan：仅液压力作用下的斜盘力；
%tolerance：允许计算误差；
%uInitial：初始解，一般输入上一步的解；
%输出参数：
%ex1、ey1、ex2、ey2：柱塞姿态；
%Fsw：斜盘对滑靴的作用力；
%iter：迭代次数；
%error：实际迭代误差；
%h：油膜厚度矩阵；
%pOil：油膜压力场矩阵；
%tautheta、tauz：油膜剪应力矩阵；
%pContact：接触应力矩阵；
%PFriction：摩擦功率；

    %解非线性方程P(u)=0
    %u为迭代过程中方程的解，每一步的解为一列向量
    %f为迭代过程中当前解对应的P(u)值，每一步的值为一列向量
    %R为剩余误差，即R=0-P(u)，每一步的值为一列向量
    %Deltau为迭代步长，Deltau(i)表示u(i+1)-u(i)，每一步的值为一列向量
    %Ks为切刚度矩阵，实时覆盖
    
    h0 = piston.h0;
    Fhtan = pHigh * pi * basic.rk^2 / cos(basic.beta_rad);
    %% 定义初始步
    u(:,1) = [uInitial(1:4) / h0 ; uInitial(5) / Fhtan];
    [f(:,1) , h(:,:,1) , pOil(:,:,1) , tautheta(:,:,1) , tauz(:,:,1) , pContact(:,:,1) , PFriction(1,1)] = ...
        BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
        pHigh , pLow , u(1) * h0 , u(2) * h0 , u(3) * h0 , u(4) * h0 , u(5) * Fhtan , phi_rad );
    R(:,1) = -f(:,1);
    
    %% 初始化误差与迭代次数
    iter = 0;
    error = norm( R(:,1) );
    
    %% 迭代循环
    while ( error(iter+1) >= tolerance && iter <= 100 )
        iter = iter + 1;
        
        %计算切刚度矩阵
        uForKs = diag( 1e-3 * ones(5,1) ) + u(:,iter);
        for i = 1:5
            [f0(:,i) , ~ , ~ , ~ , ~ , ~ , ~] = ...
                BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
                pHigh , pLow , uForKs(1,i) * h0 , uForKs(2,i) * h0 , uForKs(3,i) * h0 , uForKs(4,i) * h0 , uForKs(5,i) * Fhtan , phi_rad );
            Ks(:,i) = ( f0(:,i) - f(:,iter) ) / 1e-3;
        end
        
        Deltau(:,iter) = Ks \ R(:,iter);                %计算迭代步长
        u(:,iter + 1) = u(:,iter) + Deltau(:,iter);     %计算下一步的解    
        [f(:,iter+1) , h(:,:,iter+1) , pOil(:,:,iter+1) , tautheta(:,:,iter+1) , tauz(:,:,iter+1) , pContact(:,:,iter+1) , PFriction(iter+1,1)] = ...
            BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
            pHigh , pLow , u(1,iter+1) * h0 , u(2,iter+1) * h0 , u(3,iter+1) * h0 , u(4,iter+1) * h0 , u(5,iter+1) * Fhtan , phi_rad );
        R(:,iter + 1) = -f(:,iter + 1);
        error(iter+1,1) = norm( R(:,iter + 1) );        %计算误差
    end
    
    %% 解赋值给输出
    [error , pos] = min(error);
    ex1 = u(1,pos) * h0;
    ey1 = u(2,pos) * h0;
    ex2 = u(3,pos) * h0;
    ey2 = u(4,pos) * h0;
    Fsw = u(5,pos) * Fhtan;
    h = h(:,:,pos);
    pOil = pOil(:,:,pos);
    tautheta = tautheta(:,:,pos);
    tauz = tauz(:,:,pos);
    pContact = pContact(:,:,pos);
    PFriction = (pos);
end


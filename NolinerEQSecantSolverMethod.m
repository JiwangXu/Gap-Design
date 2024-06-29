function [ex1 , ey1 , ex2 , ey2 , Fsw , iter , error , h , pOil , tautheta , tauz , pContact , PFriction] =...
    NolinerEQSecantSolverMethod( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
    pHigh , pLow , phi_rad , tolerance , uInitial1 , uInitial2)
%求任意位置的柱塞姿态和斜盘力
%使用割线法迭代求解，相比切线法（牛顿法），割线法不需要计算雅克比矩阵，可节约大量计算时间，但迭代步数会增加；同时割线法找初始解容易不收敛；
%输入参数：
%phi_rad：当前工作角度；
%tolerance：允许计算误差，程序中误差为合力/力矩向量的模；
%uInitial、uInitial2：前两步的解，因为问题的非线性很强，因此使用前两步的解作为初始解输入可以更快找到收敛方向；
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
    %Ks为切/割刚度矩阵，实时覆盖
    
    h0 = piston.h0;
    Fhtan = pHigh * pi * basic.rk^2 / cos(basic.beta_rad);
    %% 定义第一步(初始步)
    u(:,1) = [uInitial1(1:4) / h0 ; uInitial1(5) / Fhtan];
    [f(:,1) , h(:,:,1) , pOil(:,:,1) , tautheta(:,:,1) , tauz(:,:,1) , pContact(:,:,1) , PFriction(1,1)] = ...
        BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
        pHigh , pLow , u(1) * h0 , u(2) * h0 , u(3) * h0 , u(4) * h0 , u(5) * Fhtan , phi_rad );
    R(:,1) = -f(:,1);
    
    %% 定义第二步
    u(:,2) = [uInitial2(1:4) / h0 ; uInitial2(5) / Fhtan];
    Deltau = u(:,2) - u(:,1);
    
    %如果第二步与第一步一样，则重新定义第二步为第一步所有向量往正方向走1e-3，这样可以近似计算第一步的切刚度矩阵
    if norm(Deltau) == 0
        Deltau = 0.001 * ones(5,1);
        u(:,2) = u(:,1) + Deltau;
    end
    
    [f(:,2) , h(:,:,2) , pOil(:,:,2) , tautheta(:,:,2) , tauz(:,:,2) , pContact(:,:,2) , PFriction(2,1)] = ...
        BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
        pHigh , pLow , u(1,2) * h0 , u(2,2) * h0 , u(3,2) * h0 , u(4,2) * h0 , u(5,2) * Fhtan , phi_rad );
    R(:,2) = -f(:,2);
    
    %计算初始的割刚度矩阵
    uForKs = diag( Deltau ) + u(:,1);
    for i = 1:5
        [f0(:,i) , ~ , ~ , ~ , ~ , ~ , ~] = ...
            BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
            pHigh , pLow , uForKs(1,i) * h0 , uForKs(2,i) * h0 , uForKs(3,i) * h0 , uForKs(4,i) * h0 , uForKs(5,i) * Fhtan , phi_rad );
        Ks(:,i) = ( f0(:,i) - f(:,1) ) / Deltau(i);
    end
    
    %% 初始化误差与迭代次数
    iter = 1;
    error = norm( R(:,2) );
    
    %% 迭代循环
    while ( error(iter) >= tolerance && iter <= 100)
        iter = iter + 1;
        Deltau(:,iter) = Ks \ R(:,iter);                %计算迭代步长
        u(:,iter + 1) = u(:,iter) + Deltau(:,iter);     %计算下一步的解    
        [f(:,iter+1) , h(:,:,iter+1) , pOil(:,:,iter+1) , tautheta(:,:,iter+1) , tauz(:,:,iter+1) , pContact(:,:,iter+1) , PFriction(iter+1,1)] = ...
            BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
            pHigh , pLow , u(1,iter+1) * h0 , u(2,iter+1) * h0 , u(3,iter+1) * h0 , u(4,iter+1) * h0 , u(5,iter+1) * Fhtan , phi_rad );
        R(:,iter + 1) = -f(:,iter + 1);
        Ks = Ks + ( f(:,iter + 1) - f(:,iter)  - Ks * Deltau(:,iter) ) / norm( Deltau(:,iter) )^2 * ( Deltau(:,iter) )';    %计算下一步的割刚度矩阵
        error(iter,1) = norm( R(:,iter + 1) );          %计算误差
    end
    
    %% 解赋值给输出
    [error , pos] = min(error);
    ex1 = u(1,pos+1) * h0;
    ey1 = u(2,pos+1) * h0;
    ex2 = u(3,pos+1) * h0;
    ey2 = u(4,pos+1) * h0;
    Fsw = u(5,pos+1) * Fhtan;
    h = h(:,:,pos+1);
    pOil = pOil(:,:,pos+1);
    tautheta = tautheta(:,:,pos+1);
    tauz = tauz(:,:,pos+1);
    pContact = pContact(:,:,pos+1);
    PFriction = (pos+1);
end


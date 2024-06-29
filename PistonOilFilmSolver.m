clear all ; clc ;%#ok<*SAGROW>
%% 输入泵参数
PistonParameter;
nS = 0;                         %柱塞自旋转速，正方向与主轴转动方向一致
omigaS = nS / 60 * 2 * pi;

%% 定义计算范围（角度）
phiStart = 90;
phiEnd = 270;
phiInterval = 1;

%% 网格划分密度
m = 60;         %轴向网格数
n = 60;         %角度方向网格数

%% 定义误差
tolerance = 0.1;

i = 0;
for phi = phiStart : phiInterval : phiEnd
    i = i+1;
    phi_rad(i,1) = phi * pi/180;
    
    %% 识别求解域
    [zin(i,1) , zout(i,1)] = SolveAreaCalc(basic.RT , basic.beta_rad , piston.LpIn , piston.LpOut , piston.LbIn , piston.LbOut , phi_rad(i,1) );
    zin(i,1) = zin(i,1) - basic.RT * tan(basic.beta_rad) * sin( phi_rad(i,1) );
    zout(i,1) = zout(i,1) - basic.RT * tan(basic.beta_rad) * sin( phi_rad(i,1) );
    
    %% 划分网格
    [theta_rad(:,:,i) , z(:,:,i) , m , n] = SolveAreaMesh( zin(i,1) , zout(i,1) , 0 , 2*pi , m , n );
        
    %% 计算斜盘与滑靴速度，参考系为滑靴参考系，矢径方向为为斜盘水平方向
    [vtheta1_POxyz_s , w1_POxyz_s , vtheta2_POxyz_s , w2_POxyz_s , vtheta1_COxyz_s , w1_COxyz_s , vtheta2_COxyz_s , w2_COxyz_s] =...
     VelocityCalc(basic.rk , basic.RT , basic.beta_rad , basic.omiga , omigaS , phi_rad(i,1) );
    vtheta1_POxyz(:,:,i) = vtheta1_POxyz_s * ones(m+1 , n+1);
    w1_POxyz(:,:,i) = w1_POxyz_s * ones(m+1 , n+1);
    vtheta2_POxyz(:,:,i) = vtheta2_POxyz_s * ones(m+1 , n+1);
    w2_POxyz(:,:,i) = w2_POxyz_s * ones(m+1 , n+1);
    vtheta1_COxyz(:,:,i) = vtheta1_COxyz_s * ones(m+1 , n+1);
    w1_COxyz(:,:,i) = w1_COxyz_s * ones(m+1 , n+1);
    vtheta2_COxyz(:,:,i) = vtheta2_COxyz_s * ones(m+1 , n+1);
    w2_COxyz(:,:,i) = w2_COxyz_s * ones(m+1 , n+1);

    %油膜方向的速度
    vr1(:,:,i) = zeros(m+1 , n+1);
    vr2(:,:,i) = zeros(m+1 , n+1);
    
    %% 解非线性方程，求柱塞姿态以及相关力学参数
    if i == 1
        [ex1(i,1) , ey1(i,1) , ex2(i,1) , ey2(i,1) , Fsw(i,1) , iter(i,1) , error(i,1) , ...
            h(:,:,i) , pOil(:,:,i) , tautheta(:,:,i) , tauz(:,:,i) , pContact(:,:,i) , PFriction(i,1)] =...
            NolinerEQNewtonSolverMethod( basic , physical , piston , ...
            theta_rad(:,:,i) , z(:,:,i) , vr1(:,:,i) , vr2(:,:,i) , ...
            vtheta1_COxyz(:,:,i) , vtheta2_COxyz(:,:,i) , w1_COxyz(:,:,i) , w2_COxyz(:,:,i) ,...
            basic.pr_high , basic.pr_low , phi_rad(i,1) , tolerance , zeros(5,1) );
    elseif i == 2
        [ex1(i,1) , ey1(i,1) , ex2(i,1) , ey2(i,1) , Fsw(i,1) , iter(i,1) , error(i,1) , ...
            h(:,:,i) , pOil(:,:,i) , tautheta(:,:,i) , tauz(:,:,i) , pContact(:,:,i) , PFriction(i,1)] =...
            NolinerEQNewtonSolverMethod( basic , physical , piston , ...
            theta_rad(:,:,i) , z(:,:,i) , vr1(:,:,i) , vr2(:,:,i) , ...
            vtheta1_COxyz(:,:,i) , vtheta2_COxyz(:,:,i) , w1_COxyz(:,:,i) , w2_COxyz(:,:,i) ,...
            basic.pr_high , basic.pr_low , phi_rad(i,1) , tolerance , [ex1(i-1,1) ; ey1(i-1,1) ; ex2(i-1,1) ; ey2(i-1,1) ; Fsw(i-1,1)] );
    else
%         [ex1(i,1) , ey1(i,1) , ex2(i,1) , ey2(i,1) , Fsw(i,1) , iter(i,1) , error(i,1) , ...
%             h(:,:,i) , pOil(:,:,i) , tautheta(:,:,i) , tauz(:,:,i) , pContact(:,:,i) , PFriction(i,1)] =...
%             NolinerEQNewtonSolverMethod( basic , physical , piston , ...
%             theta_rad(:,:,i) , z(:,:,i) , vr1(:,:,i) , vr2(:,:,i) , ...
%             vtheta1_COxyz(:,:,i) , vtheta2_COxyz(:,:,i) , w1_COxyz(:,:,i) , w2_COxyz(:,:,i) ,...
%             basic.pr_high , basic.pr_low , phi_rad(i,1) , tolerance , [ex1(i-1,1) ; ey1(i-1,1) ; ex2(i-1,1) ; ey2(i-1,1) ; Fsw(i-1,1)] );
        [ex1(i,1) , ey1(i,1) , ex2(i,1) , ey2(i,1) , Fsw(i,1) , iter(i,1) , error(i,1) , ...
            h(:,:,i) , pOil(:,:,i) , tautheta(:,:,i) , tauz(:,:,i) , pContact(:,:,i) , PFriction(i,1)] =...
            NolinerEQSecantSolverMethod( basic , physical , piston , ...
            theta_rad(:,:,i) , z(:,:,i) , vr1(:,:,i) , vr2(:,:,i) , ...
            vtheta1_COxyz(:,:,i) , vtheta2_COxyz(:,:,i) , w1_COxyz(:,:,i) , w2_COxyz(:,:,i) ,...
            basic.pr_high , basic.pr_low , phi_rad(i,1) , tolerance , [ex1(i-2,1) ; ey1(i-2,1) ; ex2(i-2,1) ; ey2(i-2,1) ; Fsw(i-2,1)] , ...
            [ex1(i-1,1) ; ey1(i-1,1) ; ex2(i-1,1) ; ey2(i-1,1) ; Fsw(i-1,1)] ); 
    end
    
end    
    
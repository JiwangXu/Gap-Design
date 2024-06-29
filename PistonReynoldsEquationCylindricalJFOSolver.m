function pOil = PistonReynoldsEquationCylindricalJFOSolver(pr_in , pr_out , miu0 ,...
    rk , theta_rad , z , h , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2)
%求解柱塞副压力场（雷诺方程）
%在柱坐标系下求解，其中坐标原点为柱塞球头的中心，z方向从柱塞尾部指向头部，无惯性力
%求解方法为有限差分法，采用五点式中心差分格式求解
%输入输出参数信息如下：
%pr_in、pr_out：柱塞内外的压强；
%pho0、miu0：油液密度和动力粘度[kg/m/s]；
%DeltaZ、DeltaTheta_rad：轴向及角度方向的步长；
%z：柱坐标系中网格节点的轴向坐标矩阵，第一个指标为轴向，第二个指标为角度方向；
%h：柱坐标系中网格节点的油膜厚度矩阵，第一个指标为轴向，第二个指标为角度方向；
%vr1、vr2：柱坐标系中网格节点的柱塞和缸体径向（油膜方向）速度矩阵，1为柱塞，2为缸体；
%vtheta1、vtheta2：柱坐标系中网格节点的柱塞和缸体切向速度矩阵，1为柱塞，2为缸体；
%w1、w2：柱塞和缸体轴向的运动速度矩阵，1为柱塞，2为缸体；
%pOil：输出油膜压力场矩阵。
    %% 初始化
    beta = 1.7 * 1e9;
    m = size(z,1) - 1;      %轴向网格数
    n = size(z,2) - 1;      %角度方向网格数
    Ag = zeros(m+1,n); Aconstant = zeros(m+1,n); A = zeros(m+1,n);       %p(i,j-1)系数
    Bg = zeros(m+1,n); Bconstant = zeros(m+1,n); B = zeros(m+1,n);       %p(i,j)系数
    Cg = zeros(m+1,n); Cconstant = zeros(m+1,n); C = zeros(m+1,n);       %p(i,j+1)系数
    Dg = zeros(m+1,n); Dconstant = zeros(m+1,n); D = zeros(m+1,n);       %p(i-1,j)系数
    Eg = zeros(m+1,n); Econstant = zeros(m+1,n); E = zeros(m+1,n);       %p(i+1,j)系数
    F = zeros(m+1,n);
    g = ones(m+1,n+1);
    
    %% 计算关于p(i,j)的方程的各项系数
    for i = 2:m
        for j = 2:n
            %轴向步长
            DeltaZ1 = z(i,j) - z(i-1,j);
            DeltaZ2 = z(i+1,j) - z(i,j);
            %角度方向步长
            DeltaTheta_rad1 = theta_rad(i,j) - theta_rad(i,j-1);
            DeltaTheta_rad2 = theta_rad(i,j+1) - theta_rad(i,j);
            %p(i,j-1)系数
            Ag(i,j) = 1 / rk^2 / 12 / miu0 * beta * ( ( h(i,j+1)^3 - h(i,j-1)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 * (-1) +...
                h(i,j)^3 * 2 / ( DeltaTheta_rad1 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
            Aconstant(i,j) = 1/2/rk * h(i,j-1) * ( vtheta1(i,j-1) + vtheta2(i,j-1) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 );
            A(i,j) = Aconstant(i,j) + g(i,j) * Ag(i,j); 
            %p(i,j)系数
            Bg(i,j) = 1 / 12 / miu0 * beta * h(i,j)^3 * (-2) / ( DeltaZ1 * DeltaZ2 ) +...
                1 / rk^2 / 12 / miu0 * beta * h(i,j)^3 * (-2) / ( DeltaTheta_rad1 * DeltaTheta_rad2 );
            Bconstant(i,j) = w1(i,j) * ( h(i+1,j) - h(i-1,j) ) / ( DeltaZ1 + DeltaZ2 ) +...
                vtheta2(i,j) / rk * ( h(i,j+1) - h(i,j-1) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) +...
                ( vr1(i,j) - vr2(i,j) );
            B(i,j) = Bconstant(i,j) + g(i,j) * Bg(i,j);
            %p(i,j+1)系数
            Cg(i,j) = 1 / rk^2 / 12 / miu0 * beta * ( ( h(i,j+1)^3 - h(i,j-1)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 +...
                h(i,j)^3 * 2 / ( DeltaTheta_rad2 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
            Cconstant(i,j) = -1/2/rk * h(i,j+1) * ( vtheta1(i,j+1) + vtheta2(i,j+1) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 );
            C(i,j) = Cconstant(i,j) + g(i,j) * Cg(i,j); 
            %p(i-1,j)系数
            Dg(i,j) = 1 / 12 / miu0 * beta * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 * (-1) +...
                h(i,j)^3 * 2 / ( DeltaZ1 * ( DeltaZ1 + DeltaZ2 ) ) );
            Dconstant(i,j) = 1/2 * h(i-1,j) * ( w1(i-1,j) + w2(i-1,j) ) / ( DeltaZ1 + DeltaZ2 );
            D(i,j) = Dconstant(i,j) + g(i,j) * Dg(i,j);
            %p(i+1,j)系数
            Eg(i,j) = 1 / 12 / miu0 * beta * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 +...
                h(i,j)^3 * 2 / ( DeltaZ2 * ( DeltaZ1 + DeltaZ2 ) ) );
            Econstant(i,j) = -1/2 * h(i+1,j) * ( w1(i+1,j) + w2(i+1,j) ) / ( DeltaZ1 + DeltaZ2 );
            E(i,j) = Econstant(i,j) + g(i,j) * Eg(i,j);
            %常数项系数
            F(i,j) = ( h(i+1,j) * (w1(i+1,j) + w2(i+1,j)) - h(i-1,j) * (w1(i-1,j) + w2(i-1,j)) ) / 2 / ( DeltaZ1 + DeltaZ2 ) +...
                1 / rk * ( h(i,j+1) * (vtheta1(i,j+1) + vtheta2(i,j+1)) - h(i,j-1) * (vtheta1(i,j-1) + vtheta2(i,j-1)) ) / 2 / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) -...
                w1(i,j) * ( h(i+1,j) - h(i-1,j) ) / ( DeltaZ1 + DeltaZ2 ) -...
                vtheta2(i,j) / rk * ( h(i,j+1) - h(i,j-1) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) +...
                ( vr2(i,j) - vr1(i,j) );
        end
        
        %角度方向的边缘网格节点特殊处理（j=1的情况）
        j = 1;
        DeltaZ1 = z(i,j) - z(i-1,j);
        DeltaZ2 = z(i+1,j) - z(i,j);
        DeltaTheta_rad1 = theta_rad(i,j) - theta_rad(i,n) + 2*pi;
        DeltaTheta_rad2 = theta_rad(i,j+1) - theta_rad(i,j);
        %p(i,j-1)系数
        Ag(i,j) = 1 / rk^2 / 12 / miu0 * beta * ( ( h(i,j+1)^3 - h(i,n)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 * (-1) +...
            h(i,j)^3 * 2 / ( DeltaTheta_rad1 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
        Aconstant(i,j) = 1/2/rk * h(i,n) * ( vtheta1(i,n) + vtheta2(i,n) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 );
        A(i,j) = Aconstant(i,j) + g(i,j) * Ag(i,j); 
        %p(i,j)系数
        Bg(i,j) = 1 / 12 / miu0 * beta * h(i,j)^3 * (-2) / ( DeltaZ1 * DeltaZ2 ) +...
            1 / rk^2 / 12 / miu0 * beta * h(i,j)^3 * (-2) / ( DeltaTheta_rad1 * DeltaTheta_rad2 );
        Bconstant(i,j) = w1(i,j) * ( h(i+1,j) - h(i-1,j) ) / ( DeltaZ1 + DeltaZ2 ) +...
            vtheta2(i,j) / rk * ( h(i,j+1) - h(i,n) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) +...
            ( vr1(i,j) - vr2(i,j) );
        B(i,j) = Bconstant(i,j) + g(i,j) * Bg(i,j);
        %p(i,j+1)系数
        Cg(i,j) = 1 / rk^2 / 12 / miu0 * beta * ( ( h(i,j+1)^3 - h(i,n)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 +...
            h(i,j)^3 * 2 / ( DeltaTheta_rad2 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
        Cconstant(i,j) = -1/2/rk * h(i,j+1) * ( vtheta1(i,j+1) + vtheta2(i,j+1) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 );
        C(i,j) = Cconstant(i,j) + g(i,j) * Cg(i,j); 
        %p(i-1,j)系数
        Dg(i,j) = 1 / 12 / miu0 * beta * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 * (-1) +...
            h(i,j)^3 * 2 / ( DeltaZ1 * ( DeltaZ1 + DeltaZ2 ) ) );
        Dconstant(i,j) = 1/2 * h(i-1,j) * ( w1(i-1,j) + w2(i-1,j) ) / ( DeltaZ1 + DeltaZ2 );
        D(i,j) = Dconstant(i,j) + g(i,j) * Dg(i,j);
        %p(i+1,j)系数
        Eg(i,j) = 1 / 12 / miu0 * beta * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 +...
            h(i,j)^3 * 2 / ( DeltaZ2 * ( DeltaZ1 + DeltaZ2 ) ) );
        Econstant(i,j) = -1/2 * h(i+1,j) * ( w1(i+1,j) + w2(i+1,j) ) / ( DeltaZ1 + DeltaZ2 );
        E(i,j) = Econstant(i,j) + g(i,j) * Eg(i,j);
        %常数项系数
        F(i,j) = ( h(i+1,j) * (w1(i+1,j) + w2(i+1,j)) - h(i-1,j) * (w1(i-1,j) + w2(i-1,j)) ) / 2 / ( DeltaZ1 + DeltaZ2 ) +...
            1 / rk * ( h(i,j+1) * (vtheta1(i,j+1) + vtheta2(i,j+1)) - h(i,n) * (vtheta1(i,n) + vtheta2(i,n)) ) / 2 / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) -...
            w1(i,j) * ( h(i+1,j) - h(i-1,j) ) / ( DeltaZ1 + DeltaZ2 ) -...
            vtheta2(i,j) / rk * ( h(i,j+1) - h(i,n) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) +...
            ( vr2(i,j) - vr1(i,j) );
    end

    %% 解线性方程 
    tolerance = 1/beta;
    omiga = 1;
    error = tolerance * 1e10;
    IterationNumber = 0;
    %初始化压力矩阵
    pc = 0;
    pOil = zeros(m+1,n+1);  
    pOil(1,:) = ( (pr_in - pc) / beta ) * ones(1,n+1);
    pOil(m+1,:) = ( (pr_out - pc) / beta ) * ones(1,n+1);
    
    while (error >= tolerance * 1e1)
        IterationNumber = IterationNumber + 1;  
        pOil0 = pOil;
        for i = 2:m
            j = 1;
            pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,n) - C(i,j) * pOil(i,j+1) -...
                D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
                for j = 2:n-1
                    pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,j-1) - C(i,j) * pOil(i,j+1) -...
                        D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
                end
            j = n;
            pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,j-1) - C(i,j) * pOil(i,1) -...
                D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
        end
        error = sum( sum( abs( pOil - pOil0 ) ) );
    end
    
    while (error >= tolerance)
        IterationNumber = IterationNumber + 1;
        pOil0 = pOil;
        
        g = (pOil >= 0);
        A = Aconstant + g(:,1:end-1) .* Ag;
        B = Bconstant + g(:,1:end-1) .* Bg;
        C = Cconstant + g(:,1:end-1) .* Cg;
        D = Dconstant + g(:,1:end-1) .* Dg;
        E = Econstant + g(:,1:end-1) .* Eg;
        
        for i = 2:m
            j = 1;
            if B(i,j)
                pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,n) - C(i,j) * pOil(i,j+1) -...
                    D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
            end
                for j = 2:n-1
                    if B(i,j)
                        pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,j-1) - C(i,j) * pOil(i,j+1) -...
                            D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
                    end
                end
            j = n;
            if B(i,j)
                pOil(i,j) = (1 - omiga) * pOil0(i,j) + omiga * ( F(i,j) - A(i,j) * pOil(i,j-1) - C(i,j) * pOil(i,1) -...
                    D(i,j) * pOil(i-1,j) - E(i,j) * pOil(i+1,j) ) / B(i,j);
            end
        end
        error = sum( sum( abs( pOil - pOil0 ) ) );
    end
    
    pOil = pc + g .* beta .* pOil;
    pOil(:,end) = pOil(:,1);
end
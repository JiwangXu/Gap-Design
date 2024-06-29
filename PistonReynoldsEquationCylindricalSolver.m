function pOil = PistonReynoldsEquationCylindricalSolver(pr_in , pr_out , miu0 ,...
    rk , theta_rad , z , h , hm , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2)
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
    m = size(z,1) - 1;      %轴向网格数
    n = size(z,2) - 1;      %角度方向网格数
    A = zeros(m+1,n);       %p(i,j-1)系数
    B = zeros(m+1,n);       %p(i,j)系数
    C = zeros(m+1,n);       %p(i,j+1)系数
    D = zeros(m+1,n);       %p(i-1,j)系数
    E = zeros(m+1,n);       %p(i+1,j)系数
    F = zeros(m+1,n);       %常数项系数
    K = sparse(n*(m+1) , n*(m+1));      %线性方程的刚度矩阵
    b = zeros( n*(m+1) , 1 );           %线性变化后生成的向量
    pvector = zeros( n*(m+1) , 1 );     %线性方程的解向量
    
    h(h<hm) = hm;
    
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
            A(i,j) = 1 / rk^2 / 12 / miu0 * ( ( h(i,j+1)^3 - h(i,j-1)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 * (-1) +...
                h(i,j)^3 * 2 / ( DeltaTheta_rad1 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
            %p(i,j)系数
            B(i,j) = 1 / 12 / miu0 * h(i,j)^3 * (-2) / ( DeltaZ1 * DeltaZ2 ) +...
                1 / rk^2 / 12 / miu0 * h(i,j)^3 * (-2) / ( DeltaTheta_rad1 * DeltaTheta_rad2 );
            %p(i,j+1)系数
            C(i,j) = 1 / rk^2 / 12 / miu0 * ( ( h(i,j+1)^3 - h(i,j-1)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 +...
                h(i,j)^3 * 2 / ( DeltaTheta_rad2 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
            %p(i-1,j)系数
            D(i,j) = 1 / 12 / miu0 * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 * (-1) +...
                h(i,j)^3 * 2 / ( DeltaZ1 * ( DeltaZ1 + DeltaZ2 ) ) );
            %p(i+1,j)系数
            E(i,j) = 1 / 12 / miu0 * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 +...
                h(i,j)^3 * 2 / ( DeltaZ2 * ( DeltaZ1 + DeltaZ2 ) ) );
            %线性方程等号右侧参数
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
        A(i,j) = 1 / rk^2 / 12 / miu0 * ( ( h(i,j+1)^3 - h(i,n)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 * (-1) +...
            h(i,j)^3 * 2 / ( DeltaTheta_rad1 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
        B(i,j) = 1 / 12 / miu0 * h(i,j)^3 * (-2) / ( DeltaZ1 * DeltaZ2 ) +...
            1 / rk^2 / 12 / miu0 * h(i,j)^3 * (-2) / ( DeltaTheta_rad1 * DeltaTheta_rad2 );
        C(i,j) = 1 / rk^2 / 12 / miu0 * ( ( h(i,j+1)^3 - h(i,n)^3 ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 )^2 +...
            h(i,j)^3 * 2 / ( DeltaTheta_rad2 * ( DeltaTheta_rad1 + DeltaTheta_rad2 ) ) );
        D(i,j) = 1 / 12 / miu0 * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 * (-1) +...
            h(i,j)^3 * 2 / ( DeltaZ1 * ( DeltaZ1 + DeltaZ2 ) ) );
        E(i,j) = 1 / 12 / miu0 * ( ( h(i+1,j)^3 - h(i-1,j)^3 ) / ( DeltaZ1 + DeltaZ2 )^2 +...
            h(i,j)^3 * 2 / ( DeltaZ2 * ( DeltaZ1 + DeltaZ2 ) ) );
        F(i,j) = ( h(i+1,j) * (w1(i+1,j) + w2(i+1,j)) - h(i-1,j) * (w1(i-1,j) + w2(i-1,j)) ) / 2 / ( DeltaZ1 + DeltaZ2 ) +...
            1 / rk * ( h(i,j+1) * (vtheta1(i,j+1) + vtheta2(i,j+1)) - h(i,n) * (vtheta1(i,n) + vtheta2(i,n)) ) / 2 / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) -...
            w1(i,j) * ( h(i+1,j) - h(i-1,j) ) / ( DeltaZ1 + DeltaZ2 ) -...
            vtheta2(i,j) / rk * ( h(i,j+1) - h(i,n) ) / ( DeltaTheta_rad1 + DeltaTheta_rad2 ) +...
            ( vr2(i,j) - vr1(i,j) );
    end
    
    %% Gauss-Seidel超松弛迭代求解
%     %初始化压力矩阵
%     pOil = zeros(m+1,n+1);  
%     pOil(1,:) = pr_in * ones(1,n+1);
%     pOil(m+1,:) = pr_out * ones(1,n+1);
%     [pOil , IterationNumber] = SolverGaussSeidel(A , B , C , D , E , F , pOil , 0.1 , 1.4);
%     pOil(:,end) = pOil(:,1);
    
    %% 线性方程转化为矩阵与向量的乘积格式
    for i = 2:m
        for j = 1:n
            K( (i-1)*n+j , (i-1)*n+j-1 ) = A(i,j);
            K( (i-1)*n+j , (i-1)*n+j ) = B(i,j);
            K( (i-1)*n+j , (i-1)*n+j+1 ) = C(i,j);
            K( (i-1)*n+j , (i-2)*n+j ) = D(i,j);
            K( (i-1)*n+j , i*n+j ) = E(i,j);
            b( (i-1)*n+j ) = F(i,j);
        end
        %处理边界节点
        j = 1;
        K( (i-1)*n+j , (i-1)*n+j-1 ) = 0;
        K( (i-1)*n+j , (i-1)*n+j-1+n ) = A(i,j);
        j = n;
        K( (i-1)*n+j , (i-1)*n+j+1 ) = 0;
        K( (i-1)*n+j , (i-1)*n+j+1-n ) = C(i,j);
    end

    %% 解线性方程
%     %罚方法求解
%     PenaltyCoefficient = max(max(K)) * 1e7;     %定义罚参数 
%     %已知压力的项增加罚系数
%     for i = 1:n
%         K(i,i) = PenaltyCoefficient + K(i,i);
%         b(i) = PenaltyCoefficient * pr_in;
%     end
%     for i = n*m+1:n*(m+1)
%         K(i,i) = PenaltyCoefficient + K(i,i);
%         b(i) = PenaltyCoefficient * pr_out;
%     end
%     pvector = K\b;      %解方程

    %分解法求解
    pvector(1:n) = pr_in;
    pvector( m*n+1 : n*(m+1) ) = pr_out;
    f = ( b(n+1 : m*n) - K(n+1:m*n , 1:n) * pvector(1:n) - K( n+1:m*n , m*n+1:n*(m+1) ) * pvector(m*n+1:n*(m+1)) );
    pvector(n+1 : m*n) = K(n+1:m*n , n+1:m*n) \ f;
    
    %% 解向量转化为矩阵格式
    pOil = zeros(m+1 , n);
    for i = 1:(m+1)
        pOil(i,:) = pvector( n*(i-1)+1 : n*i )';
    end
    pOil(:,end+1) = pOil(:,1);     

    pOil( pOil < 0 ) = 0;
end
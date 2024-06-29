function [Q , Q1 , Q2] = PistonOilFilmLeakgeCalc( miu0 , rk , theta_rad , z , h , w1 , w2 , pOil)
%柱塞副泄漏计算
%输入参数：
%miu0：油液动力粘度[kg/m/s]；
%rk：柱塞半径；
%theta_rad、z、h：柱坐标系中油膜的位置坐标与油膜厚度矩阵
%w1、w2：柱坐标系中内外两个面的z速度，1为内侧，2为外侧；
%pOil：油膜压力场矩阵；
%输出参数：
%Q：柱塞副总泄漏，方向为缸孔内侧往外侧的泄露量；
%Q1：压力梯度项的泄漏；
%Q2：轴向速度项产生的泄漏；
%由于输入流速场为离散形式，这里使用高斯积分法计算
%此处使用的是一维单节点高斯积分法
    %% 初始化
    m = size(z,1) - 1;      %广义x方向节点数
    n = size(z,2) - 1;      %广义y方向节点数
    W = [1;1];              %高斯积分的权系数
    xi =  [-1/sqrt(3) ; 1/sqrt(3)];    %高斯节点的xi坐标
    N = @(xi) [1-xi , 1+xi] / 2;      %形函数
    NElement = [N(xi(1)) ; N(xi(2))];    %形函数在高斯节点的值
    
    %% 单元流量计算函数
    %输入参数：
    %thetavector、hvrvector：节点的theta坐标与对应节点流量对角度的偏导值，单元为两节点，因此需输入两个节点的数据，以向量形式输入；
    %输出参数：单元流量
    function QElement = ElementQ(thetavector , hvrvector)
        JacobianElement = ( max(thetavector,[],1) - min(thetavector,[],1) ) / 2;    %雅克比行列式
        hvInterp =  NElement * hvrvector;       %高斯节点的插值
        QElement = sum(W .* JacobianElement .* hvInterp);       %单元输出结果计算
    end

    %% 泄漏计算
    %初始化
    Q1 = 0;
    Q2 = 0;
    %遍历每个单元计算泄漏流量
    for j = 1:n
        DeltaZ1 = z(m,j) - z(m-1,j);
        DeltaZ2 = z(m+1,j) - z(m,j);
        DeltaZ3 = z(m,j+1) - z(m-1,j+1);
        DeltaZ4 = z(m+1,j+1) - z(m,j+1);
        thetavector = [theta_rad(m+1,j) ; theta_rad(m+1,j+1)];
        hvrvector1(1,:) = -h(m+1,j)^3 * rk / 12 / miu0 * ( pOil(m-1,j) * DeltaZ2 / DeltaZ1 / (DeltaZ1 + DeltaZ2) -...
            pOil(m,j) * (DeltaZ1 + DeltaZ2) / DeltaZ1 / DeltaZ2 + pOil(m+1,j) * (DeltaZ1 + 2*DeltaZ2) / DeltaZ2 / (DeltaZ1 + DeltaZ2) );
        hvrvector1(2,:) = -h(m+1,j+1)^3 * rk / 12 / miu0 * ( pOil(m-1,j+1) * DeltaZ4 / DeltaZ3 / (DeltaZ3 + DeltaZ4) -...
            pOil(m,j+1) * (DeltaZ3 + DeltaZ4) / DeltaZ3 / DeltaZ4 + pOil(m+1,j+1) * (DeltaZ3 + 2*DeltaZ4) / DeltaZ4 / (DeltaZ3 + DeltaZ4) );
        hvrvector2(1,:) = h(m+1,j) * rk * ( w2(m+1,j) + w1(m+1,j) ) / 2;
        hvrvector2(2,:) = h(m+1,j+1) * rk * ( w2(m+1,j+1) + w1(m+1,j+1) ) / 2;
        QElement1 = ElementQ(thetavector , hvrvector1);
        Q1 = Q1 + QElement1;
        QElement2 = ElementQ(thetavector , hvrvector2);
        Q2 = Q2 + QElement2;
    end
    Q = Q1 + Q2;
end
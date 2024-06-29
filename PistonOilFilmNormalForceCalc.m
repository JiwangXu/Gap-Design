function [FOilx , FOily , MOilx , MOily] = PistonOilFilmNormalForceCalc(rk , theta_rad , z , pOil)
%计算油膜产生的法向液压力与力的作用位置
%输入参数：
%rk：柱塞半径；
%theta_rad、z、pOil：柱坐标系中油膜的位置坐标与压力分布矩阵
%输出参数：
%FOilx、FOily：产生的总液压支撑力
%MOilx、MOily：液压支撑力对原点产生的力矩
%由于输入压力场为离散形式，这里使用高斯积分法计算
%此处使用的是二维4节点高斯积分法
    %% 初始化
    m = size(z,1) - 1;      %广义x方向节点数
    n = size(z,2) - 1;      %广义y方向节点数
    W = [1;1;1;1];          %高斯积分的权系数
    xi =  [-1/sqrt(3) ; 1/sqrt(3) ; 1/sqrt(3) ; -1/sqrt(3)];    %高斯节点的xi坐标
    eta = [-1/sqrt(3) ; -1/sqrt(3) ; 1/sqrt(3) ; 1/sqrt(3)];    %高斯节点的eta坐标
    N4Q = @(xi , eta) [(1-xi)*(1-eta) , (1+xi)*(1-eta) , (1+xi)*(1+eta) , (1-xi)*(1+eta)] / 4;      %形函数
    N4QElement = [N4Q( xi(1),eta(1) ); N4Q( xi(2),eta(2) ) ; N4Q( xi(3),eta(3) ) ; N4Q( xi(4),eta(4) )];    %形函数在高斯节点的值
    GN4Q = @(xi , eta) [-(1-eta) , (1-eta) , (1+eta) , -(1+eta);
                        -(1-xi)  , -(1+xi) , (1+xi)  ,  (1-xi)] / 4;        %形函数的梯度
    
    %% 单元合力与合力矩计算函数
    %输入参数：
    %rk：柱塞半径
    %thetavector、zvector、pvector：节点的theta、z坐标与对应的压力值，单元为四节点，因此需输入四个节点的数据，以向量形式输入；
    %输出参数：
    %FxElement、FyElement、MxElement、MyElement：单元的合力、单元合力与单元合力对原点的力矩。
    function [FxElement , FyElement , MxElement , MyElement] = ElementForce(rk , thetavector , zvector , pvector)
        %雅克比行列式
        JacobianElement = [det( GN4Q(xi(1),eta(1) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(2),eta(2) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(3),eta(3) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(4),eta(4) ) * [zvector thetavector] )];
        %高斯节点的插值
        FxInterp =  rk * N4QElement * (-pvector .* cos(thetavector) );
        FyInterp =  rk * N4QElement * (-pvector .* sin(thetavector) );
        MxInterp =  rk * N4QElement * ( pvector .* sin(thetavector) .* zvector);
        MyInterp =  rk * N4QElement * (-pvector .* cos(thetavector) .* zvector);
        %单元输出结果计算
        FxElement = sum(W .* JacobianElement .* FxInterp);
        FyElement = sum(W .* JacobianElement .* FyInterp);
        MxElement = sum(W .* JacobianElement .* MxInterp);
        MyElement = sum(W .* JacobianElement .* MyInterp);
    end

    %% 油膜产生的液压力与力作用点计算
    %初始化，计算中心油池的液压力、力与力作用点坐标的乘积
    FOilx = 0; FOily = 0; MOilx = 0; MOily = 0;
    %遍历每个单元计算总液压力与力作用点
    for i = 1:m
        for j = 1:n
            %调取单元各个节点的坐标与函数值
            thetavector = [theta_rad(i,j+1) ; theta_rad(i,j) ; theta_rad(i+1,j) ; theta_rad(i+1,j+1)];
            zvector = [z(i,j+1) ; z(i,j) ; z(i+1,j) ; z(i+1,j+1)];
            pvector = [pOil(i,j+1) ; pOil(i,j) ; pOil(i+1,j) ; pOil(i+1,j+1)];
            %计算单元合力与合力矩
            [FxElement , FyElement , MxElement , MyElement] = ElementForce(rk , thetavector , zvector , pvector);
            %求和
            FOilx = FOilx + FxElement;
            FOily = FOily + FyElement;
            MOilx = MOilx + MxElement;
            MOily = MOily + MyElement;
        end
    end
end
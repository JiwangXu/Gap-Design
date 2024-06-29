function [pContact , Fcx , Fcy , Mcx , Mcy , Fcfx , Fcfy , Fcfz , Mcfx , Mcfy , Mcfz , PFriction] =...
    PistonContactForceCalc(rk , theta_rad , z , h , hm , vtheta , w , Episton , Ecylinder , Hcylinder , f)
%计算柱塞与缸体接触产生的作用力与作用力矩
%输入参数：
%rk：柱塞半径；
%theta_rad、z、h：柱坐标系中油膜的位置坐标与厚度矩阵；
%hm：发生接触的最大油膜厚度；
%Episton、Ecylinder：柱塞和缸体的材料弹性模量；
%Hcylinder：缸体柱塞腔的径向厚度;
%f：摩擦系数
%输出参数：
%pContact：接触应力分布矩阵；
%Fcx、Fcy：产生的总法向接触力；
%Mcx、Mcy：法向接触力对原点的作用力矩；
%Fcfx、Fcfy、Fcfz：产生的接触摩擦力；
%Mcfx、Mcfy、Mcfz：接触摩擦力对原点的作用力矩；
%PFriction：接触摩擦力的摩擦功率；
%由于输入油膜厚度场为离散形式，这里使用高斯积分法计算
%此处使用的是二维4节点高斯积分法
    %% 初始化
    Deltah = hm - h;
    Deltah(Deltah<0) = 0;
    pContact = Deltah / (2 * rk / Episton + Hcylinder / Ecylinder);
        
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
    %thetavector、zvector、Deltahvector：节点的theta、z坐标与对应的油膜厚度，单元为四节点，因此需输入四个节点的数据，以向量形式输入；
    %输出参数：
    %FcxElement、FcyElement、McxElement、McyElement：单元的合法向接触力与其对原点的力矩；
    %FcfxElement、FcfyElement、FcfzElement、McfxElement、McfyElement、McfzElement：单元的合接触摩擦力与其对原点的力矩；
    %PElement：单元的合摩擦功率。
    function [FcxElement , FcyElement , McxElement , McyElement ,...
            FcfxElement , FcfyElement , FcfzElement , McfxElement , McfyElement , McfzElement , PElement] =...
            ElementForce(thetavector , zvector , Deltahvector , vthetavector , wvector )
        %雅克比行列式
        JacobianElement = [det( GN4Q(xi(1),eta(1) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(2),eta(2) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(3),eta(3) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(4),eta(4) ) * [zvector thetavector] )];
        %高斯节点的插值
        FcInterp = Deltahvector / (2 * rk / Episton + Hcylinder / Ecylinder);
        FcxInterp =  -FcInterp .* cos(thetavector);
        FcyInterp =  -FcInterp .* sin(thetavector);
        McxInterp =   FcInterp .* sin(thetavector) .* zvector;
        McyInterp =  -FcInterp .* cos(thetavector) .* zvector;
        FcfxInterp =  vthetavector ./ sqrt(vthetavector.^2 + wvector.^2) .* FcInterp * f .* sin(thetavector);
        FcfyInterp = -vthetavector ./ sqrt(vthetavector.^2 + wvector.^2) .* FcInterp * f .* cos(thetavector);
        FcfzInterp = -wvector ./ sqrt(vthetavector.^2 + wvector.^2) .* FcInterp * f;
        McfxInterp =  rk * sin(thetavector) .* FcfzInterp - zvector .* FcfyInterp;
        McfyInterp = -rk * cos(thetavector) .* FcfzInterp + zvector .* FcfxInterp;
        McfzInterp =  rk * cos(thetavector) .* FcfyInterp - rk * sin(thetavector) .* FcfxInterp;
        PInterp = -FcInterp * f .* sqrt( vthetavector.^2 + wvector.^2 );
        
        %单元输出结果计算
        FcxElement = sum(W .* JacobianElement .* (N4QElement * FcxInterp) );
        FcyElement = sum(W .* JacobianElement .* (N4QElement * FcyInterp) );
        McxElement = sum(W .* JacobianElement .* (N4QElement * McxInterp) );
        McyElement = sum(W .* JacobianElement .* (N4QElement * McyInterp) );
        FcfxElement = sum(W .* JacobianElement .* (N4QElement * FcfxInterp) );
        FcfyElement = sum(W .* JacobianElement .* (N4QElement * FcfyInterp) );
        FcfzElement = sum(W .* JacobianElement .* (N4QElement * FcfzInterp) );
        McfxElement = sum(W .* JacobianElement .* (N4QElement * McfxInterp) );
        McfyElement = sum(W .* JacobianElement .* (N4QElement * McfyInterp) );
        McfzElement = sum(W .* JacobianElement .* (N4QElement * McfzInterp) );
        PElement = sum(W .* JacobianElement .* (N4QElement * PInterp) );
    end

    %% 柱塞与缸体接触产生的法向作用力与作用力矩计算
    %初始化
    Fcx = 0; Fcy = 0; Mcx = 0; Mcy = 0;
    Fcfx = 0; Fcfy = 0; Fcfz = 0; Mcfx = 0; Mcfy = 0; Mcfz = 0; 
    PFriction = 0;
    %遍历每个单元计算合力与合力矩
    for i = 1:m
        for j = 1:n
            %调取单元各个节点的坐标与函数值
            thetavector = [theta_rad(i,j+1) ; theta_rad(i,j) ; theta_rad(i+1,j) ; theta_rad(i+1,j+1)];
            zvector = [z(i,j+1) ; z(i,j) ; z(i+1,j) ; z(i+1,j+1)];
            Deltahvector = [Deltah(i,j+1) ; Deltah(i,j) ; Deltah(i+1,j) ; Deltah(i+1,j+1)];
            vthetavector = [vtheta(i,j+1) ; vtheta(i,j) ; vtheta(i+1,j) ; vtheta(i+1,j+1)];
            wvector = [w(i,j+1) ; w(i,j) ; w(i+1,j) ; w(i+1,j+1)];
            %计算单元合力与合力矩
            [FcxElement , FcyElement , McxElement , McyElement ,...
                FcfxElement , FcfyElement , FcfzElement , McfxElement , McfyElement , McfzElement , PElement] =...
                ElementForce(thetavector , zvector , Deltahvector , vthetavector , wvector );
            %求和
            Fcx = Fcx + FcxElement;
            Fcy = Fcy + FcyElement;
            Mcx = Mcx + McxElement;
            Mcy = Mcy + McyElement;
            Fcfx = Fcfx + FcfxElement;
            Fcfy = Fcfy + FcfyElement;
            Fcfz = Fcfz + FcfzElement;
            Mcfx = Mcfx + McfxElement;
            Mcfy = Mcfy + McfyElement;
            Mcfz = Mcfz + McfzElement;
            PFriction = PFriction + PElement;
        end
    end
end
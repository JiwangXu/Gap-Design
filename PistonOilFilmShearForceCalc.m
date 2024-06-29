function [tautheta , tauz , Ffx , Ffy , Ffz , Mfx , Mfy , PFriction] =...
    PistonOilFilmShearForceCalc( miu0 , rk , theta_rad , z , h , hm , vtheta1 , vtheta2 , w1 , w2 , pOil )
%计算油膜剪力场、摩擦力和摩擦功
%输入参数：
%miu0：油液动力粘度[kg/m/s]；
%rk：柱塞半径；
%theta_rad、z、h：柱坐标系中油膜的位置坐标与油膜厚度矩阵
%vtheta1、vtheta2、w1、w2：柱坐标系中内外两个面的theta和z速度，1为内侧，2为外侧；
%pOil：油膜压力场矩阵；
%输出参数：
%tautheta、tauz：柱坐标系中的柱塞表面剪应力场；
%Ffx、Ffy、Ffz：剪应力合力，即油液对柱塞的摩擦力；
%Mfx、Mfy、Mfz：油液对柱塞的摩擦力对原点的矩；
%PFriction：摩擦力产生的摩擦功[W]；
%由于输入参数为离散形式，这里使用高斯积分法计算
%此处使用的是二维4节点高斯积分法
    %% 初始化
    m = size(z,1) - 1;              %轴向网格数
    n = size(z,2) - 1;              %角度方向网格数
    tautheta = zeros(m+1,n+1);      %theta方向剪应力
    tauz = zeros(m+1,n+1);          %z方向剪应力
    
    h(h<hm) = hm;
    
    %% 切应力计算
    %z方向剪应力计算（中心差分）
    for i = 2:m
        DeltaZ = z(i+1,:) - z(i-1,:);
        tauz(i,:) = -h(i,:) / 2 .* ( pOil(i+1,:) - pOil(i-1,:) ) ./ DeltaZ + miu0 * ( w2(i,:) - w1(i,:) ) ./ h(i,:);
    end
    %边缘点进行特殊处理
    i = 2;
    DeltaZ1 = z(i,:) - z(i-1,:);
    DeltaZ2 = z(i+1,:) - z(i,:);
    tauz(i-1,:) = -h(i-1,:) / 2 .* ( pOil(i-1,:) .* (-2*DeltaZ1 - DeltaZ2) ./ DeltaZ1 ./ (DeltaZ1 + DeltaZ2) +...
        pOil(i,:) .* (DeltaZ1 + DeltaZ2) ./ DeltaZ1 ./ DeltaZ2 + pOil(i+1,:) .* (-DeltaZ1) ./ DeltaZ2 ./ (DeltaZ1 + DeltaZ2) ) +...
        miu0 * ( w2(i-1,:) - w1(i-1,:) ) ./ h(i-1,:);
    
    i = m;
    DeltaZ1 = z(i,:) - z(i-1,:);
    DeltaZ2 = z(i+1,:) - z(i,:);
    tauz(i+1,:) = -h(i+1,:) / 2 .* ( pOil(i-1,:) .* DeltaZ2 ./ DeltaZ1 ./ (DeltaZ1 + DeltaZ2) +...
        pOil(i,:) .* (-DeltaZ1 - DeltaZ2) ./ DeltaZ1 ./ DeltaZ2 + pOil(i+1,:) .* (DeltaZ1 + 2*DeltaZ2) ./ DeltaZ2 ./ (DeltaZ1 + DeltaZ2) ) +...
        miu0 * ( w2(i+1,:) - w1(i+1,:) ) ./ h(i+1,:);
    
    %theta方向剪应力计算
    for j = 2:n
        DeltaTheta_rad = theta_rad(:,j+1) - theta_rad(:,j-1);
        tautheta(:,j) = -h(:,j) ./ rk / 2 .* ( pOil(:,j+1) - pOil(:,j-1) ) ./ DeltaTheta_rad +...
            miu0 * ( vtheta2(:,j) - vtheta1(:,j) ) ./ h(:,j);
    end
    %边缘点进行特殊处理
    j = 1;
    DeltaTheta_rad = theta_rad(:,j+1) - theta_rad(:,j-1+n) + 2*pi;
    tautheta(:,j) = -h(:,j) ./ rk / 2 .* ( pOil(:,j+1) - pOil(:,j-1+n) ) ./ DeltaTheta_rad +...
            miu0 * ( vtheta2(:,j) - vtheta1(:,j) ) ./ h(:,j);
    tautheta(:,n+1) = tautheta(:,1);
    
    %% 初始化积分参数
    W = [1;1;1;1];          %高斯积分的权系数
    xi =  [-1/sqrt(3) ; 1/sqrt(3) ; 1/sqrt(3) ; -1/sqrt(3)];    %高斯节点的xi坐标
    eta = [-1/sqrt(3) ; -1/sqrt(3) ; 1/sqrt(3) ; 1/sqrt(3)];    %高斯节点的eta坐标
    N4Q = @(xi , eta) [(1-xi)*(1-eta) , (1+xi)*(1-eta) , (1+xi)*(1+eta) , (1-xi)*(1+eta)] / 4;      %形函数
    N4QElement = [N4Q( xi(1),eta(1) ); N4Q( xi(2),eta(2) ) ; N4Q( xi(3),eta(3) ) ; N4Q( xi(4),eta(4) )];    %形函数在高斯节点的值
    GN4Q = @(xi , eta) [-(1-eta) , (1-eta) , (1+eta) , -(1+eta);
                        -(1-xi)  , -(1+xi) , (1+xi)  ,  (1-xi)] / 4;        %形函数的梯度

    %% 单元合力与合力作用点计算函数
    %输入参数：
    %rk：柱塞半径；
    %thetavector、zvector、tauthetavector、tauzvector：节点的theta和z坐标与对应的切应力值，单元为四节点，因此需输入四个节点的数据，以向量形式输入；
    %vthetavector、wvector：节点的theta和z速度，单元为四节点，因此需输入四个节点的数据，以向量形式输入；
    %输出参数：
    %FfxElement、FfyElement、FfzElement：单元的合力；
    %MfxElement、MfyElement、MfzElement：单元对原点产生的合力；
    %PElement：单元的摩擦功。
    function [FfxElement , FfyElement , FfzElement , MfxElement , MfyElement , PElement] =...
            ElementForce(rk , thetavector , zvector , tauthetavector , tauzvector , vthetavector , wvector)
        %雅克比行列式
        JacobianElement = [det( GN4Q(xi(1),eta(1) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(2),eta(2) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(3),eta(3) ) * [zvector thetavector] ) ;
                           det( GN4Q(xi(4),eta(4) ) * [zvector thetavector] )];
        %高斯节点的插值
        FfxInterp = N4QElement * (-tauthetavector .* sin(thetavector) * rk);
        FfyInterp = N4QElement * ( tauthetavector .* cos(thetavector) * rk);
        FfzInterp = N4QElement * ( tauzvector * rk);
        MfxInterp = N4QElement * ( tauzvector .* sin(thetavector) * rk^2 - tauthetavector .* zvector .* cos(thetavector) * rk);
        MfyInterp = N4QElement * (-tauzvector .* cos(thetavector) * rk^2 - tauthetavector .* zvector .* sin(thetavector) * rk);
        PInterp = N4QElement * (tauzvector .* wvector + tauthetavector .* vthetavector) * rk;
        %单元输出结果计算
        FfxElement = sum(W .* JacobianElement .* FfxInterp);
        FfyElement = sum(W .* JacobianElement .* FfyInterp);
        FfzElement = sum(W .* JacobianElement .* FfzInterp);
        MfxElement = sum(W .* JacobianElement .* MfxInterp);
        MfyElement = sum(W .* JacobianElement .* MfyInterp);
        PElement = sum(W .* JacobianElement .* PInterp);
    end

    %% 油膜产生的总剪力和剪力作用位置计算
    %初始化
    Ffx = 0;
    Ffy = 0;
    Ffz = 0;
    Mfx = 0;
    Mfy = 0;
    PFriction = 0;
    %遍历每个单元计算合力与合力作用点
    for i = 1:m
        for j = 1:n
            %调取单元四个节点的坐标与函数值
            thetavector = [theta_rad(i,j+1) ; theta_rad(i,j) ; theta_rad(i+1,j) ; theta_rad(i+1,j+1)];
            zvector = [z(i,j+1) ; z(i,j) ; z(i+1,j) ; z(i+1,j+1)];
            tauthetavector = [tautheta(i,j+1) ; tautheta(i,j) ; tautheta(i+1,j) ; tautheta(i+1,j+1)];
            tauzvector = [tauz(i,j+1) ; tauz(i,j) ; tauz(i+1,j) ; tauz(i+1,j+1)];
            vthetavector = [vtheta1(i,j+1) ; vtheta1(i,j) ; vtheta1(i+1,j) ; vtheta1(i+1,j+1)];
            wvector = [w1(i,j+1) ; w1(i,j) ; w1(i+1,j) ; w1(i+1,j+1)];
            %计算单元合力、合力矩、摩擦功
            [FfxElement , FfyElement , FfzElement , MfxElement , MfyElement , PElement] =...
            ElementForce(rk , thetavector , zvector , tauthetavector , tauzvector , vthetavector , wvector);
            %求和
            Ffx = Ffx + FfxElement;
            Ffy = Ffy + FfyElement;
            Ffz = Ffz + FfzElement;
            Mfx = Mfx + MfxElement;
            Mfy = Mfy + MfyElement;
            PFriction = PFriction + PElement;
        end
    end

end
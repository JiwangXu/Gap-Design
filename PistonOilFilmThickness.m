function h = PistonOilFilmThickness(rk , h0 , ex1 , ey1 , z1 , ex2 , ey2 , z2 , z , theta_rad)
%求任意位置的油膜厚度。
%输入参数：
%r：圆柱（柱塞孔）半径
%h0：平均油膜厚度（柱塞孔单侧配合间隙）；
%ex1、ey1、z1、ex2、ey2、z2：两个已知位置的xy方向偏心距和z坐标；
%z、theta：要求油膜厚度位置的z坐标和角度；
%输出参数：油膜厚度；
%可以同时计算多个位置的油膜厚度，此时输入的z和theta为矩阵格式，且两个矩阵维度相同。
    if sum( size(z) ~= size(theta_rad) )
        h = 0;
        disp('计算油膜厚度时输入矩阵维度错误')
        return
    end
    %计算需求油膜厚度的z坐标位置处的偏心距
    ex = (z - z1) ./ (z2 - z1) .* (ex2 - ex1) + ex1;
    ey = (z - z1) ./ (z2 - z1) .* (ey2 - ey1) + ey1;
    e = sqrt(ex.^2 + ey.^2);
    
    %求偏心圆的圆心位置
    psi_rad = zeros(size(z,1) , size(z,2));
    psi_rad(ex ~= 0) = atan( ey(ex ~= 0) ./ ex(ex ~= 0) ) + pi * ( ex(ex ~= 0) < 0 );
    psi_rad(ex == 0) = pi/2;
    psi_rad( (ex == 0) & (ey < 0) ) = 3*pi/2;
    
    %用正弦定理计算需求油膜厚度的theta坐标位置处的中心距
    gamma_rad = psi_rad - theta_rad;
    l = rk .* sqrt( 1 - ( e ./ rk .* sin(gamma_rad) ).^2 ) + e .* cos(gamma_rad);
    
    h = h0 + rk - l;
end


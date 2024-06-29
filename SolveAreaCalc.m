function [zin , zout] = SolveAreaCalc(RT , beta_rad , LpIn , LpOut , LbIn , LbOut , phi_rad)
%柱塞副求解域计算
%计算柱塞与缸体存在油膜的轴向位置
%输入参数：
%RT、beta_rad：缸体分度圆半径与斜盘摆角
%LpIn、LpOut：柱塞中心至柱塞底部边缘与肩部边缘的轴向距离；
%LbIn、LbOut：缸体中心至缸孔内部台阶边缘与口部边缘的轴向距离；
%phi_rad：柱塞孔运动位置，起始位置为吸油区位移为0的位置；
%输出参数：
%zin、zout：求解域缸孔内侧与外侧的z坐标，z正方向指向主轴方向，因此zout>zin
%可同时输入多个运动位置，此时phi_rad以向量形式输入，输出亦为向量形式
    LWork = RT * sin(phi_rad) * tan(beta_rad);
    zout = ( -LpOut + LWork ) .* ( LpOut - LWork > LbOut ) +...
           ( -LbOut ) .* ( LpOut - LWork <= LbOut );
    zin = ( -LpIn + LWork ) .* ( LpIn - LWork < LbIn ) +...
          ( -LbIn ) .* ( LpIn - LWork >= LbIn );   
end


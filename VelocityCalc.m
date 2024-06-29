function [vtheta1_POxyz , w1_POxyz , vtheta2_POxyz , w2_POxyz ,...
          vtheta1_COxyz , w1_COxyz , vtheta2_COxyz , w2_COxyz] =...
          VelocityCalc(rk , RT , beta_rad , omiga , omigaS , phi_rad)
%计算柱塞在theta和z向的速度（不包含油膜方向速度计算）；
%输出参数：
%rk、RT、beta_rad：柱塞半径、缸体分度圆半径、斜盘摆角；
%omiga、omigaS：主轴旋转角速度、柱塞自选角速度（方向与主轴旋转方向一致，omigaS为0时柱塞做平动）；
%phi_rad：柱塞孔运动位置，起始位置为吸油区位移为0的位置；
%输出参数：
%POxyz为柱塞坐标系，非惯性系，其x方向指向全局坐标系x轴正方向；
%COxyz为缸体坐标系，非惯性系，其x方向从缸体中心指向柱塞孔中心；
%1表示柱塞，2表示缸体

    %POxsyszs系的速度，该参考系为柱塞参考系，r方向指向全局坐标系x轴正向，参考系做平动，1表示柱塞，2表示缸体
    w1_POxyz = 0;
    vtheta1_POxyz = omigaS * rk;
    w2_POxyz = -omiga * RT * cos(phi_rad) * tan(beta_rad);
    vtheta2_POxyz = omiga * rk;

    %COxsyszs系的速度，该参考系为缸体参考系，r方向从泵中心指向缸体孔中心，1表示柱塞，2表示缸体
    w2_COxyz = 0;
    vtheta2_COxyz = 0;
    w1_COxyz = omiga * RT * cos(phi_rad) * tan(beta_rad);
    vtheta1_COxyz = (-omiga + omigaS) * rk;
end


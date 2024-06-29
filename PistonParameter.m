%基础参数
basic.dk = 20 * 1e-3;                         %柱塞直径
basic.rk = basic.dk / 2;                        %柱塞半径
basic.Ak = 0.25*pi * basic.dk^2;              %柱塞面积
basic.DT = 81 * 1e-3;                         %缸体分度圆直径
basic.RT = basic.DT / 2;                      %柱塞分度圆半径
basic.beta = 16;                              %斜盘摆角
basic.beta_rad = basic.beta * pi/180;
basic.eY = 0 * 1e-3;                            %斜盘中心与柱塞中心的偏心距
basic.n = 2500;                               %转速，单位rpm
basic.omiga = 2*pi * basic.n / 60;            %角速度

basic.pr_low = (0 + 0.1013) * 1e6;                     %吸油区压力，绝对压力
basic.pr_high = (0 + 0.1013) * 1e6;                     %排油区压力，绝对压力

%物理性能参数
physical.K0 = 1700 * 1e6;                     %体积模量
physical.pho0 = 0.85 * 1e3;                   %常压下油液密度
physical.miu = 40 * 1e-3;                     %动力粘度
physical.niu0 = physical.miu / physical.pho0;     %运动粘度
physical.Episton = 2.1e11;
physical.Ecylinder = 1.16e11;

%柱塞参数
piston.h0 = 20 * 1e-6;
piston.LpOut = 10.828 * 1e-3;
piston.LpIn = 59.8 * 1e-3;
piston.LbOut = 13.95 * 1e-3;
piston.LbIn = 48.7 * 1e-3;
piston.ex1 = 17.2e-6;
piston.ey1 = 0e-6;
piston.z1 = -48.7 * 1e-3;
piston.ex2 = 19.2e-6;
piston.ey2 = 0e-6;
piston.z2 = -13.95 * 1e-3;
piston.m = 0.02;
piston.LA = 10 * 1e-3;
piston.Hcylinder = 2e-3;
piston.hm = 3*sqrt(0.1^2+0.4^2) * 1e-6;

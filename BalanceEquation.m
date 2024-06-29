function [f , h , pOil , tautheta , tauz , pContact , PFriction] = ...
    BalanceEquation( basic , physical , piston , theta_rad , z , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2 ,...
    pHigh , pLow , ex1 , ey1 , ex2 , ey2 , Fsw , phi_rad )
%柱塞平衡方程计算

    %% 计算油膜厚度
    sk = basic.RT * tan(basic.beta_rad) * sin(phi_rad);
    h = PistonOilFilmThickness(basic.rk , piston.h0 , ex1 , ey1 , piston.z1 - sk , ex2 , ey2 , piston.z2 - sk , z , theta_rad);

    %% 求解压力场
    pOil = PistonReynoldsEquationCylindricalSolver(pHigh , pLow , physical.miu ,...
        basic.rk , theta_rad , z , h , piston.hm , vr1 , vr2 , vtheta1 , vtheta2 , w1 , w2);
    
    %% 求油膜产生的力与力矩
    [FnxOil , FnyOil , MnxOil , MnyOil] = PistonOilFilmNormalForceCalc(basic.rk , theta_rad , z , pOil);
    
    %% 求接触力与力矩
    [pContact , Fcx , Fcy , Mcx , Mcy , Fcfx , Fcfy , Fcfz , Mcfx , Mcfy , Mcfz , PFriction1] =...
        PistonContactForceCalc(basic.rk , theta_rad , z , h , piston.hm ,...
        vtheta1-vtheta2 , w1-w2 , physical.Episton , physical.Ecylinder , piston.Hcylinder , 0.1);
    
    %% 求油膜产生的摩擦力与力矩
    [tautheta , tauz , Ffx , Ffy , Ffz , Mfx , Mfy , PFriction2] =...
        PistonOilFilmShearForceCalc( physical.miu , basic.rk , theta_rad , z , h , piston.hm , ...
        vtheta1 , vtheta2 , w1 , w2 , pOil );
    
    PFriction = PFriction1 + PFriction2;
    
    %% 平衡方程
    Fx = 0 + ( Fsw * sin(basic.beta_rad) * sin(phi_rad) ) + ( piston.m * basic.omiga^2 * basic.RT ) + FnxOil + Fcx + Fcfx + Ffx;
    Fy = 0 + ( Fsw * sin(basic.beta_rad) * cos(phi_rad) ) + 0 + FnyOil + Fcy + Fcfy + Ffy;
    Fz = pi * basic.rk^2 * basic.pr_high + ( -Fsw * cos(basic.beta_rad) ) + ( piston.m * basic.omiga^2 * basic.RT * tan(basic.beta_rad) * sin(phi_rad) ) + 0 + 0 + Fcfz + Ffz;
    Mx = 0 + 0 + 0 + MnxOil + Mcx + Mcfx + Mfx;
    My = 0 + 0 + piston.m * basic.omiga^2 * basic.RT  * (-piston.LA) + MnyOil + Mcy + Mcfy + Mfy;
    f = [Fx ; Fy ; Fz ; Mx ; My];
end


function [theta_rad , z , m , n] = SolveAreaMesh( z1 , z2 , theta1_rad , theta2_rad , m , n , zFine1 , zFine2)
%划分柱塞求解域网格，两个网格方向分别为z和theta方向，划分出的网格为矩形网格；
%z方向可以在两端细化网格，细化倍率调节为5；
%输入参数：
%z1、z2、theta1_rad、theta2_rad：z与theta方向划分网格的端点坐标；
%m、n：z与theta方向网格划分数目（不加密的情况下）；
%zFine1、zFine2：z向网格加密时，两端加密的截止位置；
%输出参数：
%theta_rad、z：theta和z网格坐标，以矩阵格式输出；
%m、n：z和theta方向最终的实际网格数目；
%不需要在z向细化网格时最后两个参数可不输入。
    DeltaZ = (z2 - z1) / m;             %z向网格步长
    DeltaTheta_rad = 2*pi / n;          %theta向网格步长
    MeshTheta_rad = theta1_rad : DeltaTheta_rad : theta2_rad;       %theta向网格
    %根据两个端部是否细化划分z向网格
    switch nargin
        case 6
            MeshZ = z1 : DeltaZ : z2;
        case 8
            FineCoeffient = 4;          %网格细化倍率
            MeshZ =[ z1 : DeltaZ / FineCoeffient : zFine1 ,...
                zFine1 + DeltaZ : DeltaZ : zFine2 - DeltaZ ,...
                zFine2 : DeltaZ / FineCoeffient : z2 ];
    end
    [theta_rad , z] = meshgrid( MeshTheta_rad , MeshZ );        %生成网格
    %计算网格划分后的网格数目
    m = size(z,1)-1;
    n = size(z,2)-1;
end


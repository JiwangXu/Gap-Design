function [b , IterationNumber] = SolverGaussSeidel(A , B , C , D , E , F , b0 , tolerance , omiga)
%使用Gauss-Seidel超松弛法迭代法解线性方程
%输入参数：
%A~F：线性方程的系数，分别对应(i,j-1)、(i,j+1)、(i-1,j)、(i+1,j)位置的系数；
%b0：初始解，必须定义，因为在边界位置缺少方程，初始解定义了边界条件；
%tolerance：允许误差，解矩阵中所有元素与上一次迭代元素的差的绝对值的和；
%omiga：超松弛因子，影响收敛性与求解速度，取值范围为1~2。
    m = size(A,1) - 1;
    n = size(A,2);
    error = tolerance;      %初始化误差
    IterationNumber = 0;    %初始化迭代次数
    b = b0;                 %初始化压力矩阵
    %迭代求解
    while (error >= tolerance)
        IterationNumber = IterationNumber + 1;
        b0 = b;
        for i = 2:m
            j = 1;
            if b(i,j) || IterationNumber<2
                b(i,j) = (1 - omiga) * b(i,j) + omiga * ( F(i,j) - A(i,j) * b(i,n) - C(i,j) * b(i,j+1) -...
                        D(i,j) * b(i-1,j) - E(i,j) * b(i+1,j) ) / B(i,j);
            end
            for j = 2:n-1
                if b(i,j) || IterationNumber>2
                    b(i,j) = (1 - omiga) * b(i,j) + omiga * ( F(i,j) - A(i,j) * b(i,j-1) - C(i,j) * b(i,j+1) -...
                        D(i,j) * b(i-1,j) - E(i,j) * b(i+1,j) ) / B(i,j);
                end
            end
        j = n;
        if b(i,j) || IterationNumber>2
            b(i,j) = (1 - omiga) * b(i,j) + omiga * ( F(i,j) - A(i,j) * b(i,j-1) - C(i,j) * b(i,1) -...
                    D(i,j) * b(i-1,j) - E(i,j) * b(i+1,j) ) / B(i,j);
        end
        end
        b(b<0) = 0;
        error = sum( sum( abs( b(b>0) - b0(b>0) ) ) );        %计算误差
    end
end


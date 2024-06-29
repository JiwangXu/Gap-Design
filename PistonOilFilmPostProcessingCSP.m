%以柱塞为参考系进行后处理
%% 绘制油膜厚度变化的图像
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(1);
set(figure(1) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%往5取整计算z和h的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
hmin = floor( min(min(min(h*1e6))) / 5 ) * 5;
hmax = ceil( max(max(max(h*1e6))) / 5) * 5; 
for i = 1:loops
    hh = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , h(:,:,i)*1e6 );
    grid on
    set(hh,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([hmin , hmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜厚度 [μm]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面油膜厚度场 [μm]','FontSize',15,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([hmin , hmax])
    colorbar
    colormap("jet")
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-油膜厚度.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-油膜厚度.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end

%% 绘制油膜厚度变化的图像（平面投影视图）
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(2);
set(figure(2) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%往5取整计算z和h的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
hmin = floor( min(min(min(h*1e6))) / 5 ) * 5;
hmax = ceil( max(max(max(h*1e6))) / 5) * 5; 
for i = 1:loops
    hh = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , h(:,:,i)*1e6 );
    grid on
    set(hh,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([hmin , hmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜厚度 [μm]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面油膜厚度场 [μm]','FontSize',20,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([hmin , hmax])
    colorbar
    colormap("jet")
    view(180,-90)
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-油膜厚度-平面投影.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-油膜厚度-平面投影.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end

%% 绘制油膜压力场变化的图像
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(3);
set(figure(3) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%计算z和pOil的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
pmin = floor( min(min(min(pOil/1e6))) );
pmax = ceil( max(max(max(pOil/1e6))) ); 
for i = 1:loops
    hp = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , pOil(:,:,i)/1e6 );
    grid on
    set(hp,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([pmin , pmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜压力 [MPa]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面油膜压力场 [MPa]','FontSize',20,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([pmin , pmax])
    colorbar
    colormap("jet")
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-油膜压力.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-油膜压力.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end

%% 绘制油膜压力场变化的图像（平面投影视图）
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(4);
set(figure(4) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%计算z和pOil的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
pmin = floor( min(min(min(pOil/1e6))) );
pmax = ceil( max(max(max(pOil/1e6))) ); 
for i = 1:loops
    hp = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , pOil(:,:,i)/1e6 );
    grid on
    set(hp,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([pmin , pmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜压力 [MPa]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面油膜压力场 [MPa]','FontSize',20,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([pmin , pmax])
    colorbar
    colormap("jet")
    view(180,-90)
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-油膜压力-平面投影.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-油膜压力-平面投影.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end

%% 绘制接触应力场变化的图像
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(5);
set(figure(5) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%计算z和pOil的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
pCmin = floor( min(min(min(pContact/1e6))) );
pCmax = ceil( max(max(max(pContact/1e6))) ); 
for i = 1:loops
    hpc = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , pContact(:,:,i)/1e6 );
    grid on
    set(hpc,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([pCmin , pCmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜压力 [MPa]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面接触应力场 [MPa]','FontSize',20,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([pCmin , pCmax])
    colorbar
    colormap("jet")
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-接触应力.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-接触应力.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end

%% 绘制接触应力场变化的图像（平面投影视图）
loops = size(phi_rad,1);    %提取需要生成的帧数目
figure(6);
set(figure(6) , 'Color' , 'White' , 'Position' , [460,140,1000,800]);   %设置图片背景颜色和图窗位置大小
%计算z和pOil的上下限
zmin = floor( min(min(min(z*1000))) / 5 ) * 5;
zmax = ceil( max(max(max(z*1000))) / 5 ) * 5;
pCmin = floor( min(min(min(pContact/1e6))) );
pCmax = ceil( max(max(max(pContact/1e6))) ); 
for i = 1:loops
    hpc = surf( z(:,:,i)*1000 , theta_rad(:,:,i)*180/pi , pContact(:,:,i)/1e6 );
    grid on
    set(hpc,'edgecolor','none')
    xlim([zmin , zmax])
    ylim([0 360])
    zlim([pCmin , pCmax])
    xlabel('z方向坐标 [mm]','FontSize',20,'FontWeight','bold')
    ylabel('θ方向坐标 [°]','FontSize',20,'FontWeight','bold')
    zlabel('油膜压力 [MPa]','FontSize',20,'FontWeight','bold')
    title('柱塞坐标系下柱塞表面接触应力场 [MPa]','FontSize',15,'FontWeight','bold','Color','k')
    set(gca , 'FontSize' , 20)
    caxis([pCmin , pCmax])
    colorbar
    colormap("jet")
    view(180,-90)
    drawnow;            %刷新当前绘图窗口
    M = getframe(gcf);  %获取当前帧的图片数据
    Im = frame2im(M);   %返回与电影帧相关的图片数据
    [A , map] = rgb2ind(Im,256);    %将RGB图片转化为索引图片
    %将当前帧图片写入gif图中
    if i == 1
        imwrite(A , map , 'result/柱塞坐标系-接触应力-平面投影.gif' , 'gif' , 'Loopcount' , Inf , 'DelayTime' , 0.04);
    else
        imwrite(A , map , 'result/柱塞坐标系-接触应力-平面投影.gif' , 'gif' , 'WriteMode' , 'append' , 'DelayTime' , 0.04);
    end
end
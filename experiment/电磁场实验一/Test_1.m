% 用超松弛迭代法求二维静电场域的电位分布
hx=25;hy=17;                    %设置网格节点数
v1=ones(hy,hx);                 %设置行列二维数组
v1(1:6,17:25)=zeros(6,9);
m=24;n=16;                      %横纵向网格数

% 边界的 Dirichlet 边界条件值
v1(1,1:16)=ones(1,16)*100; 
v1(2:7,16)=ones(6,1)*50;
v1(7,17:24)=ones(1,8)*50;
v1(2:17,1)=0;
v1(17,:)=0;

% 计算松弛因子
t1=(cos(pi/m)+cos(pi/n))/2;
w=2/(1+sqrt(1-t1*t1));

% 超松弛迭代法
v2=v1;maxt=1;t=0;                       % 初始化
k=0;
while(maxt>1e-6)                        % 由 v1 迭代，算出 v2，迭代精度为 0.000001
    k=k+1                               % 计算迭代次数
    maxt=0;
    for i=2:7                           % 从 2 到 7 行循环
        for j=2:15                      % 从 2 到 15 列循环
            v2(i,j)=v1(i,j)+(v1(i,j+1)+v1(i+1,j)+v2(i-1,j)+v2(i,j-1)-4*v1(i,j))*w/4;%拉普拉斯方程差分式
            t=abs(v2(i,j)-v1(i,j));
            if(t>maxt) maxt=t;
                end
        end
    end
    for i=8:(hy-1)                       % 从 8 到 hy-1 行循环
        for j=2:(hx-1)                   % 从 2 到 hx-1 列循环
            v2(i,j)=v1(i,j)+(v1(i,j+1)+v1(i+1,j)+v2(i-1,j)+v2(i,j-1)-4*v1(i,j))*w/4;%拉普拉斯方程差分式
            t=abs(v2(i,j)-v1(i,j));
            if(t>maxt) maxt=t;
                end
        end
    end
    v2(7:hy-1,hx)=v2(7:hy-1,hx-1);
    v1=v2
end

v1=v2(hy:-1:1,:)
subplot(1,2,1),mesh(v1)                                     % 画三维曲面图
 axis([0,25,0,17,0,100])
 subplot(1,2,2),contour(v1,50)                              % 画等电位线图

hold on
x=1:1:hx;y=1:1:hy;
[xx,yy]=meshgrid(x,y);                                      % 形成栅格
[Gx,Gy]=gradient(v1,0.6,0.6);                               % 计算梯度
quiver(xx,yy,Gx,Gy,'r')                                     % 根据梯度数据画箭头
axis([-1.5,hx+2.5,-2,20])                                   % 设置坐标边框
plot([1,1,hx,hx,1],[1,hy,hy,1,1],'k')                       % 画导体边框
text(hx/2-0.5,hy+0.4,'100V','fontsize',11);                 % 下标注
text(hx/2,0.3,'0V','fontsize',11);                          % 上标注
text(-0.3,hy/2,'0V','fontsize',11);                         % 左标注
text(hx+0.1,hy/2,'\partial\phi/\partialn=0','fontsize',11); % 右标注
text(hx/2+5,hy/2+6,'50V','fontsize',11);                    % 下标注
hold off
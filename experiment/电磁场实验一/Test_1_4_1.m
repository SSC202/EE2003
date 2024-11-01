% 220330124-电气1班-舒晟超-电磁场实验一-题目4

% 用超松弛迭代法求二维静电场域的电位分布
hx=25;hy=17;                    % 设置网格节点数
v1=ones(hy,hx);                 % 设置行列二维数组
v1(12:17,17:25)=zeros(6,9);     % 12到17行的17到25列电位为0
m=24;n=16;                      % 横纵向网格数

% 边界的 Dirichlet 边界条件值
v1(1,1:25)=0;                   % 第一行电位为0
v1(1:17,1)=0;                   % 第一列电位为0
v1(11,16:25)=ones(1,10)*100;    % 横向100V边界    
v1(11:17,16)=ones(7,1)*100;     % 纵向100V边界

% 计算加速收敛因子
t1=(cos(pi/m)+cos(pi/n))/2;
w=2/(1+sqrt(1-t1*t1));

% 超松弛迭代法
v2=v1;maxt=1;t=0;                       % 初始化
k=0;
while(maxt>1e-6)                        % 由 v1 迭代，算出 v2，迭代精度为 0.000001
    k=k+1                               % 计算迭代次数
    maxt=0;
    for i=2:10                          % 从 2 到 10 行循环
        for j=2:24                      % 从 2 到 24 列循环
            v2(i,j)=v1(i,j)+(v1(i,j+1)+v1(i+1,j)+v2(i-1,j)+v2(i,j-1)-4*v1(i,j))*w/4;%拉普拉斯方程差分式
            t=abs(v2(i,j)-v1(i,j));
            if(t>maxt) maxt=t;
                end
        end
    end
    for i=11:16                          % 从 11 到 16 行循环
        for j=2:15                       % 从 2 到 15 列循环
            v2(i,j)=v1(i,j)+(v1(i,j+1)+v1(i+1,j)+v2(i-1,j)+v2(i,j-1)-4*v1(i,j))*w/4;%拉普拉斯方程差分式
            t=abs(v2(i,j)-v1(i,j));
            if(t>maxt) maxt=t;
                end
        end
    end
    v2(17,2:15)=v2(16,2:15);             % Neumann 条件处理
    v2(2:10,25)=v2(2:10,24);             % Neumann 条件处理
    v1=v2                                % 迭代一次
end

v1=v2(hy:-1:1,:)
subplot(1,2,1),mesh(v1)                                     % 画三维曲面图
axis([0,25,0,17,0,100])
subplot(1,2,2),contour(v1,50)                               % 画等电位线图

hold on
x=1:1:hx;y=1:1:hy;
[xx,yy]=meshgrid(x,y);                                      % 形成栅格
[Gx,Gy]=gradient(v1,0.6,0.6);                               % 计算梯度
quiver(xx,yy,Gx,Gy,'r')                                     % 根据梯度数据画箭头
axis([-1.5,hx+2.5,-2,20])                                   % 设置坐标边框
plot([1,1,hx,hx,1],[1,hy,hy,1,1],'k')                       % 画导体边框
text(12.5,17.6,'0V','fontsize',11);   
text(-0.8,9,'0V','fontsize',11);  
text(17,3,'100V','fontsize',11);                          
text(19,5.5,'100V','fontsize',11);                         
text(25,12,'\partial\phi/\partialn=0','fontsize',11); 
text(7,0.5,'\partial\phi/\partialn=0','fontsize',11);                   
hold off
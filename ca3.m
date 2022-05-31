% 在Matlab中模拟病毒的传播模型
% 数值的建立
M=200;              % 尺寸大小
alpha=0.06;         % 健康者转化为潜伏者概率
beta=0.5;          % 潜伏者转化为传染者概率
gama=0.05;          % 康复概率
kersi=0.00052;        %死亡概率
city_r=20;

% 建立网格
%初始化人群，-1表示该位置没人，0表示健康人
y=rand([M M]);
for i=1:M
    for j=1:M
        if y(i,j)>0.2
            y(i,j)=0;
        else 
            y(i,j)=-1;
        end
    end
end


%初始化城市
%city_1
for i=1:M
    for j=1:M
        dxx = i-M/3;
        dyy = j-M/3;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<city_r && y(i,j)==-1 && rand <0.5)
            y(i,j)=0;
        end
    end
end
%city_2
for i=1:M
    for j=1:M
        dxx = i-M/1.5;
        dyy = j-M/1.5;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<city_r && y(i,j)==-1 && rand <0.5)
            y(i,j)=0;
        end
    end
end      
%city3
for i=1:M
    for j=1:M
        dxx = i-M/1.5;
        dyy = j-M/3;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<city_r && y(i,j)==-1 && rand <0.5)
            y(i,j)=0;
        end
    end
end   
%city_4
for i=1:M
    for j=1:M
        dxx = i-M/3;
        dyy = j-M/1.5;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<city_r && y(i,j)==-1 && rand <0.5)
            y(i,j)=0;
        end
    end
end    
 %y = zeros(M, M);    % 网格x元素的意思是: -1白色空位,0健康人绿色，1潜伏者黄色，2感染者红色；3康复者蓝色；4死亡者黑色
% 设置初始网格x，在网格的中心有一圈感染者，半径为10个细胞
%感染源1
for i=1:M
    for j=1:M
        dxx = i-M/2;
        dyy = j-M/2;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<2 )
            y(i,j)=2;
        end
    end
end
%感染源2
for i=1:M
    for j=1:M
        dxx = i-M/1.5;
        dyy = j-M/1.5;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<2 )
            y(i,j)=2;
        end
    end
end
y_c=reshape(y,1,40000);
n=[hist(y_c,unique(y_c)),0,0,0];
n_s=n(2);
n_e=n(4);
n_i=n(3);
n_r=n(5);
n_d=n(6);
% 隔离方块
% y(100:101,100:200)=-1;
% y(100:200,100:101)=-1;

%圆形市区
%A=reshape(round(rand([500 500])*-1),1,62500);
%B=reshape(round(rand([250 250])*500),1,62500);

% 定义社区，也就是最近的8个邻居
lingju = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0];
% 创造一个新的窗口
figure
hold on
% 主循环，迭代时间变量t
for t=1:1100
    % 遍历网格 x 中的所有单元格，对于索引 i 1从n 和 j 从1到 n

    for i=1:M
        for j=1:M
            % 在邻居之间来回走动，传播疾病
            for k=1:8
                i2 = i+lingju(k, 1);
                j2 = j+lingju(k, 2);
                % 检查细胞是否在网格边界内
                if ( i2>=1 && j2>=1 && i2<=M && j2<=M )
                    %如果细胞处于易感状态和邻近细胞
                    % 被感染的传播感染的概率是beta
                    if ( y(i,j)==0 && y(i2, j2)==2 )
                        if ( rand<alpha )
                            y(i,j) = 1;
                            continue;
                        end
                    end
                end
             end
            if (y(i,j)==1 && rand< beta)
                y(i,j)=2;
                continue;
            end
            % 如果被感染的人能够以伽马的概率从疾病中康复
            if ( y(i,j)==2 && rand<gama )
                y(i,j) = 3;
                continue;
            elseif ( y(i,j)==2 && rand<kersi)
                y(i,j) = 4;
                continue;
            end
        end
    end

    n_s=[n_s sum(y(:)==0)];
    n_e=[n_e sum(y(:)==1)];
    n_i=[n_i sum(y(:)==2)];
    n_r=[n_r sum(y(:)==3)];
    n_d=[n_d sum(y(:)==4)];
 
    % 动态模拟
    clf
    imagesc(y)                   % 网格展示
    title(['时刻：',num2str(t)]);
    pause(0.01)                         % 暂停0.01s
    colormap([1 1 1;0 1 0;1 1 0;1 0 0;0 0 1;0 0 0]);   % 定义-1,0、1、2,3,4分别对应的颜色
    fmat(:,t)=getframe;
    % 如果没有更多的感染者，就停止模拟
   %  if ( sum(y==2)==0 )
      %   break;
    % end

end

figure('color','white')
e=plot(1:length(n_s),n_s,'g');%每天健康人数变化
hold on
f=plot(1:length(n_e),n_e,'y');%每天潜伏者人数变化
hold on
g=plot(1:length(n_i),n_i,'r');%每天感染者人数变化
hold on
h=plot(1:length(n_r),n_r,'b');%每天免疫者人数变化
hold on
m=plot(1:length(n_d),n_d,'k');%每天死亡者人数变化
a=xlabel('time');
b=ylabel('number');
c=title('time-Number');
d=legend('Number S','Number E','Number I','Number R','Number D'); 
set(a,'unit','normalized','Position',[0.5,-0.05],'fontsize',15,'fontname','Times New Roman'); 
set(e,'LineWidth',1.5);
set(f,'LineWidth',1.5);
set(g,'LineWidth',1.5);
set(h,'LineWidth',1.5);
set(m,'LineWidth',1.5);
set(b,'unit','normalized','Position',[-0.05,0.5],'fontsize',15,'fontname','Times New Roman'); 
set(d,'unit','normalized','Position',[0.7,0.5,0.1,0.1],'fontsize',10)
xlim([0 length(n_s)])
xticks(0:200:length(n_s))
yticks([0:5000:n_s(1)])
grid on;

%

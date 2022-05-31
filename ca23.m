% 在Matlab中模拟病毒的传播模型
% 数值的建立
M=500;
mm=M*M;% 尺寸大小
beta=0.5263;          % 潜伏者转化为传染者概率
alpha=0.0816;         % 健康者转化为潜伏者概率
gama=0.0086;          % 康复概率
kersi=0.00052;        %死亡概率
iter=2000;
a=1;

  % 网格x元素的意思是: -1白色空位,0健康人绿色，1潜伏者黄色，2感染者红色；3康复者蓝色；4死亡者黑色
y=rand([M M]);
for i=1:M
    for j=1:M
        if y(i,j)>0.6
            y(i,j)=0;
        else 
            y(i,j)=-1;
        end
    end
end
 
% 设置初始网格x，在网格的中心有一圈感染者，半径为10个细胞
for i=1:M
    for j=1:M
        dxx = i-M/2;
        dyy = j-M/2;
        d = sqrt(dxx*dxx+dyy*dyy);
        if ( d<5 )
            y(i,j)=2;
        end
    end
end
y_c=reshape(y,1,mm);
n=[hist(y_c,unique(y_c)),0,0,0];
n_s=n(2);
n_e=n(4);
n_i=n(3);
n_r=n(5);
n_d=n(6);

% 定义社区，也就是最近的8个邻居
lingju = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0];
% 创造一个新的窗口
figure
hold on
% 主循环，迭代时间变量t
for t=1:iter
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
                    % 被感染的传播感染的概率是alpha
                    if ( y(i,j)==0 && y(i2, j2)==2 )
                        if ( rand<alpha )
                            y(i,j) = 1;
                        end
                    end
                end
             end
            if (y(i,j)==1 && rand< beta)
                y(i,j)=2;
            end
            % 如果被感染的人能够以伽马的概率从疾病中康复
            if ( y(i,j)==2 && rand<gama )
                y(i,j) = 3;
            elseif ( y(i,j)==2 && rand<kersi)
                y(i,j) = 4;
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
    if (t==50)
        print(gcf, '-dpng','E:\DESKTOP\传染病模型\ca2.3\Image1.jpg');   
    end
    if (t==iter/2.5)
        print(gcf, '-dpng','E:\DESKTOP\传染病模型\ca2.3\Image2.jpg');    
    end
    if (t==1000)
        print(gcf, '-dpng','E:\DESKTOP\传染病模型\ca2.3\Image3.jpg');    
    end
    pause(0.01)                         % 暂停0.01s
    colormap([1 1 1;0 1 0;1 1 0;1 0 0;0 0 1;0 0 0]);    % 定义0、1、2,3,4分别对应的颜色
    % 如果没有更多的感染者，就停止模拟
    if ( sum(y==2)==0 )
        a=a+1;
        if (a==100)
            break;
        end  
    end
    
end
%
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
xlim([0 iter])
xticks(0:200:iter)
yticks([0:n_s(1)/10:n_s(1)])
grid on;
set(gca,'YTick',[0:mm/10:mm]);%设置Y轴要显示坐标刻度
ms = max(n_s);
me = max(n_e);
mi = max(n_i);
mr = max(n_r);
md = max(n_d);
txt1 = ['max(S) is ',num2str(ms)];
txt2 = ['max(E) is ',num2str(me)];
txt3 = ['max(I) is ',num2str(mi)];
txt4 = ['max(R) is ',num2str(mr)];
txt5 = ['max(D) is ',num2str(md)];
text(iter/6,n_s(1)/1.6,txt1,'color','g','FontSize',10)
text(iter/6,n_s(1)/1.6-n_s(1)/20,txt2,'color','r','FontSize',10)
text(iter/6,n_s(1)/1.6-2*n_s(1)/20,txt3,'color','m','FontSize',10)
text(iter/6,n_s(1)/1.6-3*n_s(1)/20,txt4,'color','b','FontSize',10)
text(iter/6,n_s(1)/1.6-4*n_s(1)/20,txt5,'color','k','FontSize',10)
print(gcf, '-dpng', 'E:\DESKTOP\传染病模型\ca2.3\Image4.jpg' ); 

%

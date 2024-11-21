 
 
%pop——种群数量
%dim——问题维度
%ub——变量上界，[1,dim]矩阵
%lb——变量下界，[1,dim]矩阵
%fobj——适应度函数（指针）
%MaxIter——最大迭代次数
%Best_Pos——x的最佳值
%Best_Score——最优适应度值
clc;
clear all;
close all;
pop=50;
dim=2;
ub=[10,10];
lb=[-10,-10];
MaxIter=100;
fobj=@(x)fitness(x);%设置适应度函数
[Best_Pos,Best_Score,IterCurve]=GWO(pop,dim,ub,lb,fobj,MaxIter);
%…………………………………………绘图…………………………………………
figure(1);
plot(IterCurve,'r-','linewidth',2);
grid on;
title('灰狼迭代曲线');
xlabel('迭代次数');
ylabel('适应度值');
%…………………………………… 结果显示……………………………………
disp(['求解得到的x1，x2是:',num2str(Best_Pos(1)),' ',num2str(Best_Pos(2))]);
disp(['最优解对应的函数:',num2str(Best_Score)]);
 
 
%种群初始化函数
function x=initialization(pop,ub,lb,dim)
for i=1:pop
    for j=1:dim
        x(i,j)=(ub(j)-lb(j))*rand()+lb(j);
    end
end
end
%狼群越界调整函数
function x=BoundrayCheck(x,ub,lb,dim)
for i=1:size(x,1)
    for j=1:dim
        if x(i,j)>ub(j)
            x(i,j)=ub(j);
        end
        if x(i,j)<lb(j)
            x(i,j)=lb(j);
        end
    end
end
end
 
%适应度函数，可根据自身需要调整
function [Fitness]=fitness(x)
    Fitness=sum(x.^2);
end
 
 
%…………………………………………灰狼算法主体………………………………………
function [Best_Pos,Best_Score,IterCurve]=GWO(pop,dim,ub,lb,fobj,MaxIter)
Alpha_Pos=zeros(1,dim);%初始化Alpha狼群
Alpha_Score=inf;
Beta_Pos=zeros(1,dim);%初始化Beta狼群
Beta_Score=inf;
Delta_Pos=zeros(1,dim);%初始化化Delta狼群
Delta_Score=inf;
 
x=initialization(pop,ub,lb,dim);%初始化种群
Fitness=zeros(1,pop);%初始化适应度函数
for i=1:pop
    Fitness(i)=fobj(x(i,:));
end
[SortFitness,IndexSort]=sort(Fitness);
Alpha_Pos=x(IndexSort(1),:);
Alpha_Score=SortFitness(1);
Beta_Pos=x(IndexSort(2),:);
Beta_Score=SortFitness(2);
Delta_Pos=x(IndexSort(3),:);
Delta_Score=SortFitness(3);
Group_Best_Pos=Alpha_Pos;
Group_Best_Score=Alpha_Score;
for t=1:MaxIter
    a=2-t*((2)/MaxIter);%线性调整a的值
    for i=1:pop
        for j=1:dim
            %根据Alpha狼群更新位置X1
            r1=rand;
            r2=rand;
            A1=2*a*r1-a;%计算A1
            C1=2*r2;%计算C1
            D_Alpha=abs(C1*Alpha_Pos(j)-x(i,j));%计算种群中其它狼只与Alpha狼群的距离
            X1=Alpha_Pos(j)-A1*D_Alpha;%更新X1
            
            %根据Beta狼群更新位置X2
            r1=rand;
            r2=rand;
            A2=2*a*r1-a;%计算A2
            C2=2*r2;%计算C2
            D_Beta=abs(C2*Beta_Pos(j)-x(i,j));%计算种群中其它狼只与Beta狼群的距离
            X2=Beta_Pos(j)-A2*D_Beta;%更新X2
            
             %根据Delta狼群更新位置X3
            r1=rand;
            r2=rand;
            A3=2*a*r1-a;
            C3=2*r2;
            D_Delta=abs(C3*Delta_Pos(j)-x(i,j));%计算种群中其它狼只与BDelta狼群的距离
            X3=Delta_Pos(j)-A3*D_Delta;%更新X3
            
            x(i,j)=(X1+X2+X3)/3;%更新后的狼只位置
        end
    end
  x=BoundrayCheck(x,ub,lb,dim);%狼只越界调整
  for i=1:pop
      Fitness(i)=fobj(x(i,:));
      if Fitness(i)<Alpha_Score%替换Aplha狼
          Alpha_Score=Fitness(i);
          Alpha_Pos=x(i,:);
      end
      if Fitness(i)>Alpha_Score&&Fitness(i)<Beta_Score%替换Beta狼
          Beta_Score=Fitness(i);
          Beta_Pos=x(i,:);
      end
      if Fitness(i)>Alpha_Score&&Fitness(i)>Beta_Score&&Fitness(i)<Delta_Score%替换Delta狼
          Delta_Score=Fitness(i);
          Delta_Pos=x(i,:);
      end
  end
  Group_Best_Pos=Alpha_Pos;
  Group_Best_Score=Alpha_Score;
  IterCurve(t)=Group_Best_Score;
end
  Best_Pos=Group_Best_Pos;
  Best_Score=Group_Best_Score;
end
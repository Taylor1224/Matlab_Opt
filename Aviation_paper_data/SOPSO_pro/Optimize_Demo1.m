%% I. 清空环境
clc;
clear all;

%%  
fitnessfcn = @fun;
% 变量个数
nvars=4;
% lb<= X <= ub
lb=[0,0,0,0];
ub=[];
% A*X <= b 
A = [0    0 1 1
    -1/3  0 0 0
     0 -1/2 0 0
     0    0 0 0];
 
b = [48 ; 30 ; 30 ; 0];
 
% Aeq*X = beq
Aeq=[1 1 0 0;0 0 0 0; 0 0 0 0; 0 0 0 0];
beq=[120;0;0;0];
%最优个体系数paretoFraction
%种群大小populationsize
%最大进化代数generations
%停止代数stallGenLimit
%适应度函数偏差TolFun
%函数gaplotpareto：绘制Pareto前沿 
options=gaoptimset('paretoFraction',0.3,'populationsize',200,'generations',300,'stallGenLimit',200,'TolFun',1e-10,'PlotFcns',@gaplotpareto);
 
[x,fval] = gamultiobj(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,options)
%  x是优化问题的最优解; fval是优化问题的最优解对应的目标函数值。
%  options是一个结构体，包含了一些优化选项，例如种群大小、迭代次数等。
plot(-fval(:,1),fval(:,1),'pr')
xlabel('f_1(x)')
ylabel('f_2(x)')
title('Pareto front')
grid on
%% I. ��ջ���
clc;
clear all;

%%  
fitnessfcn = @fun;
% ��������
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
%���Ÿ���ϵ��paretoFraction
%��Ⱥ��Сpopulationsize
%����������generations
%ֹͣ����stallGenLimit
%��Ӧ�Ⱥ���ƫ��TolFun
%����gaplotpareto������Paretoǰ�� 
options=gaoptimset('paretoFraction',0.3,'populationsize',200,'generations',300,'stallGenLimit',200,'TolFun',1e-10,'PlotFcns',@gaplotpareto);
 
[x,fval] = gamultiobj(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,options)
%  x���Ż���������Ž�; fval���Ż���������Ž��Ӧ��Ŀ�꺯��ֵ��
%  options��һ���ṹ�壬������һЩ�Ż�ѡ�������Ⱥ��С�����������ȡ�
plot(-fval(:,1),fval(:,1),'pr')
xlabel('f_1(x)')
ylabel('f_2(x)')
title('Pareto front')
grid on
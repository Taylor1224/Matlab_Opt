clc;clear;

Nf1 = load('D_LIFE0.5.txt');
Nf2 = load('D_LIFE0.9.txt');
Nf3 = load('D_LIFE0.95.txt');
Nf4 = load('D_LIFE0.5_optim.txt');
Nf5 = load('D_LIFE0.9_optim.txt');
Nf6 = load('D_LIFE0.95_optim.txt');

n=[100:100:5500]';

for j = 1:1:length(n)
   z1(:,j) = 1-n(j,1)./Nf1;
   a1 = find(z1(:,j)>0);
   R1(1,j)=length(a1)/length(z1(:,j));

   z2(:,j) = 1-n(j,1)./Nf2;
   a2 = find(z2(:,j)>0);
   R2(1,j)=length(a2)/length(z2(:,j));

   z3(:,j) = 1-n(j,1)./Nf3;
   a3 = find(z3(:,j)>0);
   R3(1,j)=length(a3)/length(z3(:,j));

   z4(:,j) = 1-n(j,1)./Nf4;
   a4 = find(z4(:,j)>0);
   R4(1,j)=length(a4)/length(z4(:,j));

   z5(:,j) = 1-n(j,1)./Nf5;
   a5 = find(z5(:,j)>0);
   R5(1,j)=length(a5)/length(z5(:,j));

   z6(:,j) = 1-n(j,1)./Nf6;
   a6 = find(z6(:,j)>0);
   R6(1,j)=length(a6)/length(z6(:,j));
end

% 画图：
plot(n,R1,'--');
hold on
plot(n,R2,'--');
hold on 
plot(n,R3,'--');
hold on 
plot(n,R4,'-');
hold on 
plot(n,R5,'-');
hold on 
plot(n,R6,'-');

% 标签：
xlabel('Applied cycle n/cycles');
ylabel('Reliability R');

legend ('P=0.5 before optimization', 'P=0.9 before optimization', 'P=0.95 before optimization', 'P=0.5 after optimization','P=0.9 after optimization','P=0.95 after optimization');

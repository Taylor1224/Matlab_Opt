% clc;clear;
A1=readtable('BLISK4_D123.xlsx','sheet','Fatigue','range','O2:O2001','ReadVariableNames',false);
A=table2array(A1);
B1=readtable('BLISK4_D123.xlsx','sheet','Creep','range','X2:X2001','ReadVariableNames',false);
B=table2array(B1);
for i=1:2000
    C(i)=20800/A(i);
    D(i)=20800*(19.41-2.57)/B(i);
    E(i)=C(i)+D(i);
end
dfittool(E);
% probplot('lognormal',E);
xx01=PDFda(:,1);
uu01=PDFda(:,2);
xx02=PDFda1(:,1);
uu02=PDFda1(:,2);
xx03=PDFda2(:,1);
uu03=PDFda2(:,2);
xx04=PDFda3(:,1);
uu04=PDFda3(:,2);
xx05=PDFda4(:,1);
uu05=PDFda4(:,2);
plot(xx01,uu01,'r','LineWidth',2);
hold on
plot(xx02,uu02,'b','LineWidth',2);
plot(xx03,uu03,'g','LineWidth',2);
plot(xx04,uu04,'k','LineWidth',2);
plot(xx05,uu05,'y','LineWidth',2);
set(gca,'XTick',0:0.4:1.2);
set(gca,'XLim',[0 1.2]);
set(gca,'XTickLabel',{'0','0.4','0.8','1.2'});
xlabel('\fontname{Times New Roman}\fontsize{16}Damage {\itD_A} ','color','k');
set(gca,'FontName','Times New Roman','FontSize',16);
set(gca,'YLim',[0 25]);
set(gca,'ytick',0:5:25);
set(gca,'YTickLabel',{'0','0.0025','0.0050','0.0075','0.0100','0.0125'});
ylabel('\fontname{Times New Roman}\fontsize{16}Density','fontsize',16);
legend('3050','5200','10400','15400','20800','box',' off');
x01=PDFd(:,1);
u01=PDFd(:,2);
x02=PDFd1(:,1);
u02=PDFd1(:,2);
x03=PDFd2(:,1);
u03=PDFd2(:,2);
x04=PDFd3(:,1);
u04=PDFd3(:,2);
x05=PDFd4(:,1);
u05=PDFd4(:,2);
plot(x01,u01,'r','LineWidth',2);
hold on
plot(x02,u02,'b','LineWidth',2);
plot(x03,u03,'g','LineWidth',2);
plot(x04,u04,'k','LineWidth',2);
plot(x05,u05,'y','LineWidth',2);
set(gca,'XTick',0:0.02:0.06);
set(gca,'XLim',[0 0.06]);
set(gca,'XTickLabel',{'0','0.02','0.04','0.06'});
xlabel('\fontname{Times New Roman}\fontsize{16}Damage {\itD_B} ','color','k');
set(gca,'FontName','Times New Roman','FontSize',16);
set(gca,'YLim',[0 400]);
set(gca,'ytick',0:100:400);
set(gca,'YTickLabel',{'0','0.005','0.010','0.015','0.020'});
ylabel('\fontname{Times New Roman}\fontsize{16}Density','fontsize',16);
legend('3050','5200','10400','15400','20800','box',' off');
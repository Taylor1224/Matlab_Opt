%% 四个疲劳参数的样本：
N = 100;
% 置信度为0.5时：
muX1 = [1.5943e3;-0.1058;250.6938;-1.4698];
sigmaX1 = [79.7150;5.5941e-8;20.055504;1.9847e-5];
cvX1 = sigmaX1(3)/muX1(3);
sLn1 = sqrt(log(1+cvX1^2));mLn1 = log(muX1(3))-sLn1^2/2; % 对数正态分布的均值和方差。 
y1 = lhsdesign(N,4);
x1 = [norminv(y1(:,1),muX1(1),sigmaX1(1)),norminv(y1(:,2),muX1(2),sigmaX1(2)),logninv(y1(:,3),mLn1,sLn1),norminv(y1(:,4),muX1(4),sigmaX1(4))];
% 置信度为0.9时：
muX2 = [1.4776e3;-0.1058;76.6160;-1.4698];
sigmaX2 = [73.88;5.5941e-8;6.12928;1.9847e-5];
cvX2 = sigmaX2(3)/muX2(3);
sLn2 = sqrt(log(1+cvX2^2));mLn2 = log(muX2(3))-sLn2^2/2;
y2 = lhsdesign(N,4);
x2 = [norminv(y2(:,1),muX2(1),sigmaX2(1)),norminv(y2(:,2),muX2(2),sigmaX2(2)),logninv(y2(:,3),mLn2,sLn2),norminv(y2(:,4),muX2(4),sigmaX2(4))];
% 置信度为0.95时：
muX3 = [1.4439e3;-0.1058;53.3976;-1.4698];
sigmaX3 = [72.195;5.5941e-8;4.271808;1.9847e-5];
cvX3 = sigmaX3(3)/muX3(3);
sLn3 = sqrt(log(1+cvX3^2));mLn3 = log(muX3(3))-sLn3^2/2;
y3 = lhsdesign(N,4);
x3 = [norminv(y3(:,1),muX3(1),sigmaX3(1)),norminv(y3(:,2),muX3(2),sigmaX3(2)),logninv(y3(:,3),mLn3,sLn3),norminv(y3(:,4),muX3(4),sigmaX3(4))];


%% 弹性模量E的样本：
muX0 = 161000;
sigmaX0 = 3220;
y = lhsdesign(N,1);
x = norminv(y(:,1),muX0,sigmaX0) ;
e = x;



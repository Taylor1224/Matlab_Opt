function predy = kriging(test_x)

train_x = readmatrix('BLISK_D123.csv','Range','E1:H2000');
train_y = readmatrix('BLISK_D123.csv','Range','P1:P2000');


%无需归一化，直接赋值

nva = 4;
finaltheta = ones(1,nva);

%训练
[dmodel, perf] = dacefit(train_x, train_y, @regpoly0, @corrgauss, finaltheta);
% @regpoly0，@regpoly1，@regpoly2可选，分别是常数多项式，一次多项式，二次多项式
[predy] = predictor(test_x, dmodel);%预测



end




 


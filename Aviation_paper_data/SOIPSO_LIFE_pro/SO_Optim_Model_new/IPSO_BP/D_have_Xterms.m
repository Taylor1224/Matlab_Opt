function y = D_have_Xterms(x)
    a = 5.75725e-3;
    a0 = [1.14931e-4,4.34528e-5,-1.26323e-5,8.26056e-5];
    
    B = [4.60196e-6,1.42121e-6/2,1.00806e-6/2,-7.28609e-6/2;
    1.42121e-6/2,1.20870e-7,4.48961e-7/2,-1.26447e-6/2;
    1.00806e-6/2,4.48961e-7/2,3.39699e-8,-1.94344e-6/2;
    -7.28609e-6/2,-1.26447e-6/2,-1.94344e-6/2,1.13531e-6];

    c = [4.91966e-2,0,0,0;
    0,7.57591e-2,0,0;
    0,0,2.17238e-1,0;
    0,0,0,1.00009e-1];

    d = [5.00036e1;5.00008e1;4.99648e1;5.00045e1];

    x = [x(1);x(2);x(3);x(4)];%这里的x(i)是得出的响应面模型里的X_scaled !!!
    y = a+a0*(c*x-d)+(c*x-d)'*B*(c*x-d);

    % ！！！函数名必须与文件名相同！！！
    % ！！！函数定义时，function前面不能加任何注释（%）或者（%%）！！！
end
function y = BDstress_DOF_Xterms(x)

    c = [4.91966e-2,0,0,0;
        0,7.57591e-2,0,0;
        0,0,2.17238e-1,0;
        0,0,0,1.00009e-1];

    d = [5.00036e1;5.00008e1;4.99648e1;5.00045e1];

    x = [x(1);x(2);x(3);x(4)];

    % B的目标函数：
    a1 = 7.14057e2;

    b1 = [6.49155,1.00674e1,-3.54765,1.28289e1];
    
    B1 = [6.85135e-2,-6.99269e-2/2,-1.42584e-2/2,1.62663e-1/2;
        -6.99269e-2/2,1.14680e-1,-9.25117e-4/2,-1.29160e-1/2;
        -1.42584e-2/2,-9.25117e-4/2,-2.49299e-2,-1.18437e-2/2;
        1.62663e-1/2,-1.29160e-1/2,-1.18437e-2/2,-1.68695e-2];

    y1 = a1+b1*(c*x-d)+(c*x-d)'*B1*(c*x-d);

    % D的目标函数：
    a2 = 9.42084e2;

    b2 = [1.87925e1,7.10923,-2.295,1.09977e1];
    
    B2 = [7.52309e-1,2.36426e-1/2,1.55476e-1/2,-1.2285/2;
        2.36426e-1/2,1.91950e-2,7.75868e-2/2,-2.18328e-1/2;
        1.55476e-1/2,7.75868e-2/2,4.81757e-3,-3.56018e-1/2;
        -1.2285/2,-2.18328e-1/2,-3.56018e-1/2,-9.18835e-2];

    y2 = a2+b2*(c*x-d)+(c*x-d)'*B2*(c*x-d);

    % DOF_27387(叶片径向变形)的目标函数：
    a3 = 9.90567e-1;

    b3 = [8.53132e-3,6.16370e-3,3.89536e-3,1.16272e-2];
    
    B3 = [9.60864e-5,4.62971e-6/2,6.97340e-6/2,-1.64547e-5/2;
        4.62971e-6/2,2.32914e-5,-6.11107e-7/2,-6.30436e-6/2;
        6.97340e-6/2,-6.11107e-7/2,3.62693e-5,5.49038e-5/2;
        -1.64547e-5/2,-6.30436e-6/2,5.49038e-5/2,5.17507e-5];

    y3 = a3+b3*(c*x-d)+(c*x-d)'*B3*(c*x-d);

    % DOF_88479(盘径向变形)的目标函数：
    a4 = 3.05195e-1;

    b4 = [5.72482e-3,-3.72223e-4,3.98613e-3,1.32553e-2];
    
    B4 = [5.73221e-5,1.33162e-6/2,2.37206e-6/2,4.51975e-6/2;
        1.33162e-6/2,2.50375e-7,-1.39120e-6/2,-5.38669e-6/2;
        2.37206e-6/2,-1.39120e-6/2,3.26689e-5,5.37955e-5/2;
        4.51975e-6/2,-5.38669e-6/2,5.37955e-5/2,3.33542e-5];

    y4 = a4+b4*(c*x-d)+(c*x-d)'*B4*(c*x-d);

    % 总的目标函数：
    y = [y1 y2 y3 y4]';

end
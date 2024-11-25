function y = objective_function(x)
    % Objective function (multi-objective, e.g., ZDT1)
    c = [4.91966e-2,0,0,0;
        0,7.57591e-2,0,0;
        0,0,2.17238e-1,0;
        0,0,0,1.00009e-1];

    d = [5.00036e1;5.00008e1;4.99648e1;5.00045e1];

    x = [x(1),x(2),x(3),x(4)];

    % B的LCF目标函数：
    a1 = 7.14057e2;

    b1 = [6.49155,1.00674e1,-3.54765,1.28289e1];
    
    B1 = [6.85135e-2,-6.99269e-2/2,-1.42584e-2/2,1.62663e-1/2;
        -6.99269e-2/2,1.14680e-1,-9.25117e-4/2,-1.29160e-1/2;
        -1.42584e-2/2,-9.25117e-4/2,-2.49299e-2,-1.18437e-2/2;
        1.62663e-1/2,-1.29160e-1/2,-1.18437e-2/2,-1.68695e-2];

    y1 = a1+b1*(c*x'-d)+(c*x'-d)'*B1*(c*x'-d);

    % DOF_27387(叶片径向变形)的目标函数：
    a3 = 9.90567e-1;

    b3 = [8.53132e-3,6.16370e-3,3.89536e-3,1.16272e-2];
    
    B3 = [9.60864e-5,4.62971e-6/2,6.97340e-6/2,-1.64547e-5/2;
        4.62971e-6/2,2.32914e-5,-6.11107e-7/2,-6.30436e-6/2;
        6.97340e-6/2,-6.11107e-7/2,3.62693e-5,5.49038e-5/2;
        -1.64547e-5/2,-6.30436e-6/2,5.49038e-5/2,5.17507e-5];

    y3 = a3+b3*(c*x'-d)+(c*x'-d)'*B3*(c*x'-d);
    % 总的目标函数：
    y = [y1,y3];
end
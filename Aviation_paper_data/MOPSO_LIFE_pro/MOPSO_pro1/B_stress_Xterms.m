function y = B_stress_Xterms(x)
    a = 7.14057e2;

    a0 = [6.49155,1.00674e1,-3.54765,1.28289e1];
    
    B = [6.85135e-2,-6.99269e-2/2,-1.42584e-2/2,1.62663e-1/2;
    -6.99269e-2/2,1.14680e-1,-9.25117e-4/2,-1.29160e-1/2;
    -1.42584e-2/2,-9.25117e-4/2,-2.49299e-2,-1.18437e-2/2;
    1.62663e-1/2,-1.29160e-1/2,-1.18437e-2/2,-1.68695e-2];

    c = [4.91966e-2,0,0,0;
        0,7.57591e-2,0,0;
        0,0,2.17238e-1,0;
        0,0,0,1.00009e-1];

    d = [5.00036e1;5.00008e1;4.99648e1;5.00045e1];

    x = [x(1);x(2);x(3);x(4)];

    y = a+a0*(c*x-d)+(c*x-d)'*B*(c*x-d);
end
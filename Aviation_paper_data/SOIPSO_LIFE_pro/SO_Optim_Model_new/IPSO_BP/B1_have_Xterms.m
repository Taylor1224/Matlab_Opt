function y = B1_have_Xterms(x)
    a = 4.43321e-3;
    a0 = [4.11802e-5,7.16894e-5,-1.95189e-5,8.82892e-5];
    
    B = [-4.60746e-7,3.95316e-6/2,-1.2346e-6/2,8.15143e-6/2;
    3.95316e-6/2,-2.66953e-6,1.41298e-6/2,-6.65472e-6/2;
    -1.2346e-6/2,1.41298e-6/2,6.71814e-7,-3.57939e-6/2;
    8.15143e-6/2,-6.65472e-6/2,-3.57939e-6/2,-7.99105e-6];

    c = [5.17798e-2,0,0,0;
    0,6.89986e-2,0,0;
    0,0,2.10794e-1,0;
    0,0,0,1.09157e-1];

    d = [5.26337e1;4.55018e1;4.84582e1;5.45754e1];
    
    x = [x(1);x(2);x(3);x(4)];
    y = a+a0*(c*x-d)+(c*x-d)'*B*(c*x-d);

end
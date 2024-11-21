clc;clear;

%% Latin超立方抽样求寿命：
nS = 5e2;nS1 = 50;

%% 优化前的寿命：
muX0 = [4.4269444e-3;0.71385e3;161000;0];
sigmaX0 = [9.2887787e-5;0.13671e2;3220;1];
NF1_0 = ones(50,1);
NF2_0 = ones(50,1);
NF3_0 = ones(50,1);
for k = 1:nS/nS1
    y = lhsdesign(nS1,4);
    x = [norminv(y(:,1),muX0(1),sigmaX0(1)),norminv(y(:,2),muX0(2),sigmaX0(2)),norminv(y(:,3),muX0(3),sigmaX0(3)),norminv(y(:,4),muX0(4),sigmaX0(4))];

    % 置信度为0.5：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF1_0(m,1) = Nf;
    end
    %Nf = solve((10.^(-18.956*ones(nS1,1)+0.206*x(:,4))./(9.4553*ones(nS1,1)+5.0017e-6*x(:,4))-x(:,2)./x(:,3)).*(2*Nf).^(-1*ones(nS1,1)./(9.4553*ones(nS1,1)+5.0017e-6*x(:,4))) ...
        %+ (10.^(1.542*ones(nS1,1)+0.23117*x(:,4))./(0.68038*ones(nS1,1)+9.1883e-6*x(:,4))).*(2*Nf).^(-1*ones(nS1,1)./(0.68038*ones(nS1,1)+9.1883e-6*x(:,4)))==x(:,1)/2,Nf);

    % 置信度为0.9：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF2_0(m,1) = Nf;
    end

    % 置信度为0.95：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF3_0(m,1) = Nf;
    end
end
% histfit(NF)

%% 优化后的寿命：
muX = [0.0042051;680.5919;161000;0];
sigmaX = [3.7024e-05;4.8827;3220;1];
NF1 = ones(50,1);
NF2 = ones(50,1);
NF3 = ones(50,1);
for k = 1:nS/nS1
    y = lhsdesign(nS1,4);
    x = [norminv(y(:,1),muX0(1),sigmaX0(1)),norminv(y(:,2),muX0(2),sigmaX0(2)),norminv(y(:,3),muX0(3),sigmaX0(3)),norminv(y(:,4),muX0(4),sigmaX0(4))];

    % 置信度为0.5：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF1(m,1) = Nf;
    end
    %Nf = solve((10.^(-18.956*ones(nS1,1)+0.206*x(:,4))./(9.4553*ones(nS1,1)+5.0017e-6*x(:,4))-x(:,2)./x(:,3)).*(2*Nf).^(-1*ones(nS1,1)./(9.4553*ones(nS1,1)+5.0017e-6*x(:,4))) ...
        %+ (10.^(1.542*ones(nS1,1)+0.23117*x(:,4))./(0.68038*ones(nS1,1)+9.1883e-6*x(:,4))).*(2*Nf).^(-1*ones(nS1,1)./(0.68038*ones(nS1,1)+9.1883e-6*x(:,4)))==x(:,1)/2,Nf);

    % 置信度为0.9：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF2(m,1) = Nf;
    end

    % 置信度为0.95：
    for m = 1:nS1
        syms Nf
        
        Nf = vpasolve((10^((-18.956+0.206*x(m,4))/(9.4553+5.0017e-6*x(m,4)))-x(m,2)/x(m,3))*(2*Nf)^(-1/(9.4553+5.0017e-6*x(m,4))) ...
            + 10^((1.542+0.23117*x(m,4))/(0.68038+9.1883e-6*x(m,4)))*(2*Nf)^(-1/(0.68038+9.1883e-6*x(m,4)))==x(m,1)/2,Nf);
        NF3(m,1) = Nf;
    end
end
% histfit(NF)
function error = fun1(pop, hiddennum, net, Pn_train, Tn_train)

    inputnum  = size(Pn_train, 1);   % 输入层节点数
    outputnum = size(Tn_train, 1);   % 输出层节点数
    

    w1 = pop(1 : inputnum * hiddennum);
    B1 = pop(inputnum * hiddennum + 1 : inputnum * hiddennum + hiddennum);
    w2 = pop(inputnum * hiddennum + hiddennum + 1 : ...
        inputnum * hiddennum + hiddennum + hiddennum * outputnum);
    B2 = pop(inputnum * hiddennum + hiddennum + hiddennum * outputnum + 1 : ...
        inputnum * hiddennum + hiddennum + hiddennum * outputnum + outputnum);

    net.Iw{1, 1} = reshape(w1, hiddennum, inputnum );
    net.Lw{2, 1} = reshape(w2, outputnum, hiddennum);
    net.b{1}     = reshape(B1, hiddennum, 1);
    net.b{2}     = B2';

    net = train(net, Pn_train, Tn_train);

    t_sim1 = sim(net, Pn_train);

    error = sqrt(sum((t_sim1 - Tn_train) .^ 2) ./ length(t_sim1));
end
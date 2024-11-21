x1 = linspace(-10,10,200)';
x2 = linspace(-10,10,200)';
y = exp(-x1.^2)+x2.^2;

a = [x1,x2,y];
fid = fopen('Numerical_experiment.csv','wt');%写入文件路径
[m,n]=size(a); 
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',a(i,j));
      else
        fprintf(fid,'%g\t',a(i,j));
       end
    end
end
fclose(fid);
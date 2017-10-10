function res = MyDist1(p,dim)
r = sqrt(p(1)^2+p(2)^2);
r1 = (p(1)+0.2)^2+(p(2)+0.2)^2;
res = (r-1)^2*exp(-r1/3);

function res = MyDist2(p1,dim)
% SampleDistribution: mixture of two normal distributions
res = exp( -p1*p1)+0.5*exp(-(p1+3)^2);

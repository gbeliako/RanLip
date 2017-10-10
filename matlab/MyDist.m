function res = MyDist(p,dim)
% SampleDistribution
res = exp( -(p(2)-p(1)^2)^2 -  (p(1)^2+p(2)^2)/2);

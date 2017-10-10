dim = 1;
Left = [-50];
Right = [50];

ranlip('PrepareHatFunctionAuto',dim, Left, Right, @MyDist2, 10000,8,1);
disp 'computed lipschitz constant ='  
lc = ranlip('Lipschitz')

p = ranlip('RandomVec')
disp 'Generate 50000  random vectors using ranlip("RandomVec",50000)'  

X = ranlip('RandomVec', 50000);
hist(X,40)


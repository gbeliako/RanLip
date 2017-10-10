
dim = 2;
Left = [-3,-3];
Right = [3,3];
[X,Y]=meshgrid(-3:.125:3);
disp 'Plot the distribution density function MyDist in [-3..3]x[-3..3]';
clear Z
for i=1:49 for j=1:49 Z(i,j)=MyDist([X(i,j),Y(i,j)],2); end; end;
surfc(X,Y,Z)
view([110,40])
disp 'Press any key to continue...';
pause
disp 'Next generate the hat function using ranlip("PrepareHatFunctionAuto",dim, Left, Right, @MyDist, 10,8,1);';
ranlip('PrepareHatFunctionAuto',dim, Left, Right, @MyDist, 10,8,1);
disp 'computed lipschitz constant ='  
lc = ranlip('Lipschitz')
ranlip('Seed', 10);
disp 'Next generate a couple of random vectors using ranlip("RandomVec")'  
for i = 1:2
	p = ranlip('RandomVec')
end
%errCnt  = ranlip('Count_error')
disp 'Generate 5000  random vectors using ranlip("RandomVec",5000)'  
X = ranlip('RandomVec', 5000);
figure
disp 'Scatter plot generated vectors, compare with the density function'  
plot(X(:,1),X(:,2),'.')

disp 'Save the hat function into a file...'  

ranlip('savepartition','part1.txt')

disp '   Clear the memory. We cannot generate random vectors anymore...' 
ranlip('Freemem')
p = ranlip('RandomVec')

disp 'Next load the hat function and set the distribution function' 
ranlip('loadpartition','part1.txt')
ranlip('setdistfunction',@MyDist)

disp '   Resume random vector generation' 
p = ranlip('RandomVec')

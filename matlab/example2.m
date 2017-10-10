dim = 2;
Left = [-4;-4];
Right = [4;4];
[X,Y]=meshgrid(-4:.125:4);
for i=1:65 for j=1:65 Z(i,j)=MyDist1([X(i,j),Y(i,j)],2); end; end;
surfc(X,Y,Z)
view([110,40])
disp 'Press any key to continue...';
pause
ranlip('PrepareHatFunction',dim, Left, Right, @MyDist1, 20,8,1);
p = ranlip('RandomVec')

disp 'Generate 5000  random vectors using ranlip("RandomVec",5000)'  
X = ranlip('RandomVec', 5000);
figure
disp 'Scatter plot generated vectors, compare with the density function'  
plot(X(1,:),X(2,:),'.')


function plotFM(rx,ry,A,B,N)

xx=linspace(-rx,rx,N);

yy=linspace(-ry,ry,N);

[X, Y] = meshgrid(xx,yy);

z1= A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2 ;
z2= B(1,1)*(X).^2 + 2*B(1,2)*(X).*(Y) + B(2,2)*(Y).^2;

contour(xx,yy,z1,'color','blue');
hold on;
contour(xx,yy,z2,'color','blue');
hold off;

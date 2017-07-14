function plotFM(rx,ry,A,N)

xx=linspace(-rx,rx,N);

yy=linspace(-ry,ry,N);

[X, Y] = meshgrid(xx,yy);

z =A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2;

contour(xx,yy,z,'color','red');
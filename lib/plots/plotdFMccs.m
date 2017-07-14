function plotdFMccs(rx,ry,A,B,N)

xx=linspace(-rx,rx,N);

yy=linspace(-ry,ry,N);

[X, Y] = meshgrid(xx,yy);

z1= A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2 ;
z2= B(1,1)*(X).^2 + 2*B(1,2)*(X).*(Y) + B(2,2)*(Y).^2;
%colormap(HSV);

%figure
[c, h]=contour( xx, yy, z1, 'LineWidth', 1, 'color','red');
%www=get(h,'LevelList')

%set(h,'LevelList',[10 2050 4050 6050 8050 10050]/1000)


hold on
 
[c1, h1]=contour( xx, yy, z2, 'LineWidth', 1, 'color','black');
%get(h1,'LevelList');

%set(h1,'LevelList',www)
hold off
    
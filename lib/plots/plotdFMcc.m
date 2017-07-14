function plotdFMcc(name,rx,ry,C,D,i,j,L1,L2)

[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');

A=C([i j],[i j]);
B=D([i j],[i j]);

N=100;
xx=linspace(-rx,rx,N);

yy=linspace(-ry,ry,N);

[X, Y] = meshgrid(xx,yy);

z1= A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2 ;
z2= B(1,1)*(X).^2 + 2*B(1,2)*(X).*(Y) + B(2,2)*(Y).^2;
%colormap(HSV);

[c, h]=contour( xx, yy, z1, 'LineWidth', 1, 'color','red');
get(h,'LevelList');
hold on
 
[c1, h1]=contour( xx, yy, z2, 'LineWidth', 1, 'color','black');
get(h1,'LevelList'); 
hold off

xlabel(parn(i));
ylabel(parn(j));



legend(L1, L2)
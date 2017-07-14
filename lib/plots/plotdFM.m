function plotdFM(name,rx,ry,C,D,i,j)
N=100;


[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q');

A=C([i j],[i j]);
B=D([i j],[i j]);

xx=linspace(-rx,rx,N);

yy=linspace(-ry,ry,N);

[X, Y] = meshgrid(xx,yy);

z1= A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2 ;
z2= B(1,1)*(X).^2 + 2*B(1,2)*(X).*(Y) + B(2,2)*(Y).^2;
colormap(hsv);

hSurf=pcolor(xx,yy,(z1));
shading flat;
shading interp;

hSurfAx=(gca);
cRange= caxis; % get the default color map for the data range
hold on

[C hT]= contour( xx, yy, z2, 'LineWidth', 1.5, 'color','black');
   % hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
   % set(hLines, 'black','LineWidth', 2); % and set their width.
caxis(cRange);
    
hold off
  
xlabel(parn(i));
ylabel(parn(j));
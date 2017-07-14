function plotdFMnorm(rx,ry,A,B,N)


xx=linspace(-rx/sqrt(A(1,1)),rx/sqrt(A(1,1)),N);

yy=linspace(-ry/sqrt(A(2,2)),ry/sqrt(A(2,2)),N);

[X, Y] = meshgrid(xx,yy);

z1= A(1,1)*(X).^2 + 2*A(1,2)*(X).*(Y) + A(2,2)*(Y).^2 ;
z2= B(1,1)*(X).^2 + 2*B(1,2)*(X).*(Y) + B(2,2)*(Y).^2;
colormap(HSV);

%figure
hSurf=pcolor(xx,yy,(z1));
shading flat;
shading interp;

hSurfAx=(gca);
cRange= caxis; % get the default color map for the data range
hold on

[C hT]= contour( xx, yy, z2, 'LineWidth', 0.5, 'color','black');


   % hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
   % set(hLines, 'black','LineWidth', 2); % and set their width.
caxis(cRange);
    
hold off
    
% figure
% pcolor(xx,yy,z1);
% shading flat;
% shading interp;
% 
% %contourf(xx,yy,z1,lev, 'LineWidth', 0.1);
% %set(hand, 'LineWidth', 0.1);
 %hold on;
% %figure
%contour(xx,yy,z2);
% contour(xx,yy,z2,'color','black','LineWidth', 2);

 % %set(hand, 'LineWidth', 2);
%hold off;

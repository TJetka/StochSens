function R = f(t,y,p)

R = [
     [ ((2*p(3)*y(3))/(p(4) + y(1))^2 - (2*p(3)*y(1)*y(3))/(p(4) + y(1))^3)  (0)  ((p(3)*y(1))/(p(4) + y(1))^2 - p(3)/(p(4) + y(1))) ];
     [ (0)  (0)  (0) ];
     [ (0)  (0)  (0) ];
];

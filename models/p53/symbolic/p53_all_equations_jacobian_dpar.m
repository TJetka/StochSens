function R = f(t,y,p)

R = [
     [ (1)  (-y(1))  (-(y(1)*y(3))/(p(4) + y(1)))  ((p(3)*y(1)*y(3))/(p(4) + y(1))^2)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (y(1))  (-y(2))  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (y(2))  (-y(3)) ];
     [ (1)  (y(1) - 2*y(4))  ((y(1)*y(3))/(p(4) + y(1)) - 2*y(4)*(y(3)/(p(4) + y(1)) - (y(1)*y(3))/(p(4) + y(1))^2) - (2*y(1)*y(8))/(p(4) + y(1)))  (2*y(4)*((p(3)*y(3))/(p(4) + y(1))^2 - (2*p(3)*y(1)*y(3))/(p(4) + y(1))^3) - (p(3)*y(1)*y(3))/(p(4) + y(1))^2 + (2*p(3)*y(1)*y(8))/(p(4) + y(1))^2)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (y(1) + 2*y(7))  (y(2) - 2*y(5))  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (y(2) + 2*y(9))  (y(3) - 2*y(6)) ];
     [ (0)  (-y(7))  (- y(7)*(y(3)/(p(4) + y(1)) - (y(1)*y(3))/(p(4) + y(1))^2) - (y(1)*y(9))/(p(4) + y(1)))  (y(7)*((p(3)*y(3))/(p(4) + y(1))^2 - (2*p(3)*y(1)*y(3))/(p(4) + y(1))^3) + (p(3)*y(1)*y(9))/(p(4) + y(1))^2)  (y(4))  (-y(7))  (0) ];
     [ (0)  (-y(8))  (- y(8)*(y(3)/(p(4) + y(1)) - (y(1)*y(3))/(p(4) + y(1))^2) - (y(1)*y(6))/(p(4) + y(1)))  (y(8)*((p(3)*y(3))/(p(4) + y(1))^2 - (2*p(3)*y(1)*y(3))/(p(4) + y(1))^3) + (p(3)*y(1)*y(6))/(p(4) + y(1))^2)  (0)  (y(7))  (-y(8)) ];
     [ (0)  (0)  (0)  (0)  (y(8))  (y(5) - y(2) - y(9))  (-y(9)) ];
];
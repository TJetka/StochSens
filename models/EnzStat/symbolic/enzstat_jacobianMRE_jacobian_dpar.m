function R = f(t,y,p)

R = [
     [ (-y(2))  (0)  (0)  (0)  (0) ];
     [ (-y(1))  (0)  (0)  (0)  (0) ];
     [ (0)  (1)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0) ];
     [ (-y(2))  (0)  (0)  (0)  (0) ];
     [ (-y(1))  (0)  (0)  (0)  (0) ];
     [ (0)  (1)  (1)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0) ];
     [ (y(2))  (0)  (0)  (0)  (0) ];
     [ (y(1))  (0)  (0)  (0)  (0) ];
     [ (0)  (-1)  (-1)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (1)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (-1) ];
];
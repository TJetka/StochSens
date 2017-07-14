function R = f(t,y,p)

R = [
     [ (singlegene_stimulus_stimulus(t, 1))  (0)  (-y(1))  (0) ];
     [ (0)  (y(1))  (0)  (-y(2)) ];
     [ (singlegene_stimulus_stimulus(t, 1))  (0)  (y(1) - 2*y(3))  (0) ];
     [ (0)  (y(1) + 2*y(5))  (0)  (y(2) - 2*y(4)) ];
     [ (0)  (y(3))  (-y(5))  (-y(5)) ];
];

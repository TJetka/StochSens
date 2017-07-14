function R = JAKSTATtest_rates(x, par, t, stimulus)
R =[
par(1)*stimulus(1)*x(1);
par(2)*x(2)*x(2)*stimulus(2);
par(3)*x(3)*stimulus(3);
par(4)*x(4);
];
 end

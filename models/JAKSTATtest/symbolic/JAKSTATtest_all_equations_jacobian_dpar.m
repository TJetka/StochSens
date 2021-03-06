function R = f(t,y,p)

R = [
     [ (-y(1)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1))  (-2*y(2)^2*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (y(2)^2*JAKSTATtest_stimulus(t, 2))  (-y(3)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(3)*JAKSTATtest_stimulus(t, 3))  (-y(4))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(4))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(3)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(3)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(3)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1) - 2*y(18)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1) + 2*y(35)*JAKSTATtest_stimulus(t, 1))  (4*y(2)^2*JAKSTATtest_stimulus(t, 2) - 8*y(2)*y(19)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (y(2)^2*JAKSTATtest_stimulus(t, 2) + 4*y(2)*y(51)*JAKSTATtest_stimulus(t, 2))  (y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(20)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(3)*JAKSTATtest_stimulus(t, 3) + 2*y(66)*JAKSTATtest_stimulus(t, 3))  (y(4) - 2*y(21))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(4) + 2*y(80))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1) + 2*y(47)*JAKSTATtest_stimulus(t, 1))  (0)  (4*y(3)*JAKSTATtest_stimulus(t, 3) - 4*y(76)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (4*y(3)*JAKSTATtest_stimulus(t, 3) - 4*y(77)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(3)*JAKSTATtest_stimulus(t, 3) + 2*y(78)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(18)*JAKSTATtest_stimulus(t, 1) - y(1)*JAKSTATtest_stimulus(t, 1) - y(35)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(35)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(36)*JAKSTATtest_stimulus(t, 1))  (2*y(2)*y(35)*JAKSTATtest_stimulus(t, 2))  (-y(36)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(37)*JAKSTATtest_stimulus(t, 1))  (0)  (y(36)*JAKSTATtest_stimulus(t, 3))  (-y(37))  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(38)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (y(37))  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(39)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(40)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(41)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(42)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(43)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(44)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(45)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(46)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(18)*JAKSTATtest_stimulus(t, 1) - y(1)*JAKSTATtest_stimulus(t, 1) - y(47)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(36)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(48)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(36)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(49)*JAKSTATtest_stimulus(t, 1))  (0)  (y(36)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (-y(50)*JAKSTATtest_stimulus(t, 1))  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(36)*JAKSTATtest_stimulus(t, 1))  (2*y(2)*y(19)*JAKSTATtest_stimulus(t, 2) - 2*y(2)^2*JAKSTATtest_stimulus(t, 2) - 4*y(2)*y(51)*JAKSTATtest_stimulus(t, 2))  (-y(51)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(37)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(52)*JAKSTATtest_stimulus(t, 2))  (y(51)*JAKSTATtest_stimulus(t, 3))  (-y(52))  (0)  (0)  (0)  (0)  (0) ];
     [ (y(38)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(53)*JAKSTATtest_stimulus(t, 2))  (0)  (y(52))  (0)  (0)  (0)  (0)  (0) ];
     [ (y(39)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(54)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(40)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(55)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(41)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(56)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(42)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(57)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(43)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(58)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(44)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(59)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(45)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(60)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(46)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(61)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(1)*JAKSTATtest_stimulus(t, 1) + y(35)*JAKSTATtest_stimulus(t, 1) + y(47)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(62)*JAKSTATtest_stimulus(t, 2))  (-2*y(51)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(48)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(63)*JAKSTATtest_stimulus(t, 2))  (-2*y(51)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(49)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(64)*JAKSTATtest_stimulus(t, 2))  (y(51)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(50)*JAKSTATtest_stimulus(t, 1))  (-4*y(2)*y(65)*JAKSTATtest_stimulus(t, 2))  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(52)*JAKSTATtest_stimulus(t, 2))  (y(20)*JAKSTATtest_stimulus(t, 3) - y(3)*JAKSTATtest_stimulus(t, 3) - y(66)*JAKSTATtest_stimulus(t, 3))  (-y(66))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(53)*JAKSTATtest_stimulus(t, 2))  (-y(67)*JAKSTATtest_stimulus(t, 3))  (y(66))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(54)*JAKSTATtest_stimulus(t, 2))  (-y(68)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(55)*JAKSTATtest_stimulus(t, 2))  (-y(69)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(56)*JAKSTATtest_stimulus(t, 2))  (-y(70)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(57)*JAKSTATtest_stimulus(t, 2))  (-y(71)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(58)*JAKSTATtest_stimulus(t, 2))  (-y(72)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(59)*JAKSTATtest_stimulus(t, 2))  (-y(73)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(60)*JAKSTATtest_stimulus(t, 2))  (-y(74)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(61)*JAKSTATtest_stimulus(t, 2))  (-y(75)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(36)*JAKSTATtest_stimulus(t, 1))  (2*y(2)*y(62)*JAKSTATtest_stimulus(t, 2))  (2*y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(20)*JAKSTATtest_stimulus(t, 3) - y(76)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(63)*JAKSTATtest_stimulus(t, 2))  (2*y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(20)*JAKSTATtest_stimulus(t, 3) - y(77)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(64)*JAKSTATtest_stimulus(t, 2))  (y(20)*JAKSTATtest_stimulus(t, 3) - y(3)*JAKSTATtest_stimulus(t, 3) - y(78)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (2*y(2)*y(65)*JAKSTATtest_stimulus(t, 2))  (-y(79)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(67)*JAKSTATtest_stimulus(t, 3))  (y(21) - y(4) - y(80))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(68)*JAKSTATtest_stimulus(t, 3))  (-y(81))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(69)*JAKSTATtest_stimulus(t, 3))  (-y(82))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(70)*JAKSTATtest_stimulus(t, 3))  (-y(83))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(71)*JAKSTATtest_stimulus(t, 3))  (-y(84))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(72)*JAKSTATtest_stimulus(t, 3))  (-y(85))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(73)*JAKSTATtest_stimulus(t, 3))  (-y(86))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(74)*JAKSTATtest_stimulus(t, 3))  (-y(87))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(75)*JAKSTATtest_stimulus(t, 3))  (-y(88))  (0)  (0)  (0)  (0)  (0) ];
     [ (y(37)*JAKSTATtest_stimulus(t, 1))  (0)  (y(76)*JAKSTATtest_stimulus(t, 3) - 2*y(66)*JAKSTATtest_stimulus(t, 3) - 2*y(3)*JAKSTATtest_stimulus(t, 3))  (-y(89))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(77)*JAKSTATtest_stimulus(t, 3) - 2*y(66)*JAKSTATtest_stimulus(t, 3) - 2*y(3)*JAKSTATtest_stimulus(t, 3))  (-y(90))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(3)*JAKSTATtest_stimulus(t, 3) + y(66)*JAKSTATtest_stimulus(t, 3) + y(78)*JAKSTATtest_stimulus(t, 3))  (-y(91))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(79)*JAKSTATtest_stimulus(t, 3))  (-y(92))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(81))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(82))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(83))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(84))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(85))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(86))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(87))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(88))  (0)  (0)  (0)  (0)  (0) ];
     [ (y(38)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(67)*JAKSTATtest_stimulus(t, 3))  (y(89))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(67)*JAKSTATtest_stimulus(t, 3))  (y(90))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(67)*JAKSTATtest_stimulus(t, 3))  (y(91))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (y(92))  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(39)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(68)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(68)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(68)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(40)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(69)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(69)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(69)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(41)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(70)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(70)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(70)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(42)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(71)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(71)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(71)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(43)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(72)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(72)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(72)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(44)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(73)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(73)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(73)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(45)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(74)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(74)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(74)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(46)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(75)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(75)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(75)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(48)*JAKSTATtest_stimulus(t, 1))  (0)  (4*y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(76)*JAKSTATtest_stimulus(t, 3) - 2*y(77)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(49)*JAKSTATtest_stimulus(t, 1))  (0)  (y(76)*JAKSTATtest_stimulus(t, 3) - 2*y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(78)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (y(50)*JAKSTATtest_stimulus(t, 1))  (0)  (-2*y(79)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(77)*JAKSTATtest_stimulus(t, 3) - 2*y(3)*JAKSTATtest_stimulus(t, 3) - 2*y(78)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (-2*y(79)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
     [ (0)  (0)  (y(79)*JAKSTATtest_stimulus(t, 3))  (0)  (0)  (0)  (0)  (0)  (0) ];
];

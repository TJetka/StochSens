function R = Enzyme_rates(y,par)
R=[par(4);
par(1)*y(1)*y(2);
par(2)*y(3);
par(3)*y(3);
par(5)*y(4);];
end


function R = model1_rates(y,par)
R=[par(1);
 par(2)*y(1);
 par(3)*y(1);
 par(4)*y(2);];
end


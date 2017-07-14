function R = LV_rates(y,par,t)
R=[par(1)*y(1);
par(2)*y(1)*y(2);
par(3)*y(2);];
end


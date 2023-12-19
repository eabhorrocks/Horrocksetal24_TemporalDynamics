function rho = fastPearsonCorr(x, y)

n = length(x);
x_mean = mean(x);
y_mean = mean(y);
x_std = std(x);
y_std = std(y);
rho = sum((x - x_mean) .* (y - y_mean)) / (n * x_std * y_std);

end
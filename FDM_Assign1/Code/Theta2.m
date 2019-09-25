function [y] = Theta1(x, t, rr)
    y = exp(-((x+t)/rr).^2);
end


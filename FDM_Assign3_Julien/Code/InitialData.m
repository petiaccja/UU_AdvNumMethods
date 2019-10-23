function [u1, u2] = InitialData(x, rr)
    u1 = Theta2(x, 0, rr) - Theta1(x, 0, rr);
    u2 = Theta2(x, 0, rr) + Theta1(x, 0, rr);
end
function [u1, u2] = AnalyticSolution(x, t, L, rr)
    u1 = Theta1(x, L-L*t, rr) - Theta2(x, L-L*t, rr);
    u2 = Theta1(x, L-L*t, rr) + Theta2(x, L-L*t, rr);
end
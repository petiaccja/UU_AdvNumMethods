function [u] = AnalyticSolution(x, t, c, a, eps)
     u = c - a*tanh(a*(x-c*t)/(2*eps));
end
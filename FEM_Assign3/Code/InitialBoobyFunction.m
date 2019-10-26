function z = InitialBoobyFunction(x, y)
    r_0 = 0.25;
    x1_0 = 0.3;
    x2_0 = 0;

    z = 1/2*(1-tanh(((x-x1_0)^2+(y-x2_0)^2)/(r_0^2)-1));
end
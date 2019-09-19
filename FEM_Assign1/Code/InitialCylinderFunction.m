function z = InitialCylinderFunction(x, y)
    r_0 = 0.25;
    x0 = 0.3;
    y0 = 0.0;

    if ((x-x0)^2+(y-y0)^2 <= r_0^2)
        z = 1;
    else
        z = 0;
    end
end
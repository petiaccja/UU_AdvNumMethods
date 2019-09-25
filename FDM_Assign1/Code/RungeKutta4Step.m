function [ ynext ] = RungeKutta4Step( f, tn, yn, h )
    k_1 = h*f(tn, yn);
    k_2 = h*f(tn+h/2, yn+k_1/2);
    k_3 = h*f(tn+h/2, yn+k_2/2);
    k_4 = h*f(tn+h, yn+k_3);
    
    ynext = yn + 1/6*(k_1+2*k_2+2*k_3+k_4);
end
    
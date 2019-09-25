function [ y ] = ANM_FDM_RK4( f, t_0, y_0, h, n )
%Standard Runge Kutta Method

t_old = t_0;
y_old = y_0;
y = y_0;
for i=1:n
    k_1 = h*f(t_old, y_old);
    k_2 = h*f(t_old+h/2, y_old+k_1/2);
    k_3 = h*f(t_old+h/2, y_old+k_2/2);
    k_4 = h*f(t_old+h, y_old+k_3);
    
    y_new = y_old + 1/6*(k_1+2*k_2+2*k_3+k_4);
    
    t_old = t_old + h;
    y_old = y_new;
    y = [y,y_old];
end %for, i



end


h = [1/4,1/8,1/16,1/32];

L2E = zeros(length(h), 1);

for i=1:length(h)
    %mesh
    g = @circleg;
    [p,e,t] = initmesh(g, 'hmax', h(i));
    [ U, M ] = ANM_ass1_FEM_unstable( p, e, t );
    e = U(:,1)-U(:,end);
    L2E(i) = sqrt(e'*M*e);
end %for, hmax

ls = polyfit(h', L2E, 2);
x = linspace(0,1/4,26);
y = zeros(1, 26);
for i=1:length(x)
    y(i) = ls(1)*x(i)^2+ls(2)*x(i)+ls(3);
end %i
size(x)
size(y)
plot(h, L2E, 'r', x, y, 'b', 'linewidth', 2)
legend('Error function', 'Quadratic interpol.')
xlabel('h')
ylabel('L^2-error')
print('ANM_Ass1_data1_convergence','-djpeg')


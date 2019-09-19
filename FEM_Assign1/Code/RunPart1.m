main()

% Main function
function main()
    h = [1/4, 1/8, 1/16, 1/32];

    L2E = zeros(length(h), 1);
    
    for i=1:length(h)
        geometry = @circleg;
        [p,e,t] = initmesh(geometry, 'hmax', h(i));
        
        xi = CreateInitialData(p, @InitialBoobyFunction);

        [U, M] = SolverCN(p, t, xi, 1/50, 50);

        e = U(:,1)-U(:,end);
        L2E(i) = sqrt(e'*M*e);
    end
    
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
end


% Initial datas
function z = InitialBoobyFunction(x, y)
    r_0 = 0.25;
    x1_0 = 0.3;
    x2_0 = 0;

    z = 1/2*(1-tanh(((x-x1_0)^2+(y-x2_0)^2)/(r_0^2)-1));
end

function z = InitialCylinderFunction(x, y) 
    % TODO:
    z = 0;
end
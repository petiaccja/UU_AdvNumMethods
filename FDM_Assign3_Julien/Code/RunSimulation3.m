function v = RunSimulation3(gridDim, deltaT, numberOfIterations, MakeSBPOperators)

    %initialisation
    c = 2;
    a = 1;
    eps = .00001;
    
    %space
    h = 2/(gridDim-1);
    x = linspace(-1,1,gridDim)';

    %initial data
    v = AnalyticSolution(x,0,c,a,eps);
    [HI, D1, D2, DD, M] = MakeSBPOperators(gridDim);
    Hinv = HI;
    
    gamma = 1;
    DI = -gamma*HI*DD'*DD;
    
    
    %runge kutta step
    for i=1:numberOfIterations
        %v = RungeKutta4Step(@(t,v) SAT_iterate(t,v,D1,D2,Hinv,eps), i*deltaT, v, deltaT);
        v = RungeKutta4Step(@(t,v) SAT_iterate_AD(t,v,D1,D2,DI,Hinv,eps), i*deltaT, v, deltaT);
        plot(v)
        ylim([0,5])
        drawnow
        if isnan(sum(v))
            break;
        end
    end
end
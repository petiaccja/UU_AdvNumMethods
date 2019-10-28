function v = RunSimulation3(gridDim, deltaT, numberOfIterations, MakeSBPOperators, eps)

    %initialisation
    c = 2;
    a = 1;
    
    
    %space
    h = 2/(gridDim-1);
    x = linspace(-1,1,gridDim)';

    %initial data
    v = AnalyticSolution(x,0,c,a,eps);
    plot(v)
    
    %[HI, D1, D2, DD, M, MakeD2] = MakeSBPOperators(gridDim);
    [HI, D1, D2, DD, M] = MakeSBPOperators(gridDim);
    Hinv = HI;
    
    gamma = 1;
    DI = -gamma*HI*DD'*DD;
    
    
    %runge kutta step
    %vprevaaa = v;
    for i=1:numberOfIterations
        %v = RungeKutta4Step(@(t,v) SAT_iterate(t,v,D1,D2,Hinv,eps), i*deltaT, v, deltaT);
        v = RungeKutta4Step(@(t,v) SAT_iterate_AD(t,v,D1,D2,DI,Hinv,eps), i*deltaT, v, deltaT);
        
        %Iterate = @(t,v,vprev) SAT_iterate_RV( t, v, vprev, D1, D2, MakeD2, Hinv, deltaT, h, eps);
        %v_new = RungeKutta4StepEx(Iterate, i*deltaT, v, vprevaaa, deltaT);
        %vprevaaa = v;
        %v = v_new;
        %figure(1);
        plot(v)
        ylim([0,5])
        drawnow
        if isnan(sum(v))
            break;
        end
    end
    
end
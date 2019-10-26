function Sd = ConvectionMatrixLS(p,t) 
    np=size(p,2); nt=size(t,2); 
    Sd=sparse(np,np); 
    for i=1:nt 
        loc2glb=t(1:3,i); 
        x=p(1,loc2glb); y=p(2,loc2glb); 
        [area,b,c]=HatGradients(x,y); 
        bxmid = mean(-2*pi*y);
        bymid = mean(2*pi*x);
        SdK=(bxmid*b+bymid*c)*(bxmid*b+bymid*c)'*area; 
        Sd(loc2glb,loc2glb)=Sd(loc2glb,loc2glb)+SdK; 
    end
end
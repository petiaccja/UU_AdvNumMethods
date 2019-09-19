%initial data
r_0 = 0.25;
x1_0 = 0.3;
x2_0 = 0;
T = 1;
u_0 = @(x1,x2) 1/2*(1-tanh(((x1-x1_0)^2+(x2-x2_0)^2)/(r_0^2)-1));
f = @(x, y) 0;
f_der = @(x1, x2) 2*pi*[-x2,x1]  ;
CFL = 0.5;

%hmax
hmax = 1/8;
%hmax = 1/16;

%time partition
N = 25;
time = linspace(0, 1, N+1);
delta_t = T/N;

%mesh
g = @circleg;
[p,e,t] = initmesh(g, 'hmax', hmax);


%construct the identity matrix
I = eye(length(p));
%stiffness matrix
A = AssembleA(p,e,t);
bx = @(y) 2*pi*(-y);
by = @(x) 2*pi*x;
C = ConvectionAssembler2D(p,t,bx,by);
M = MassAssembler2D(p,t);

np = size(p,2);
xi_old = zeros(np,1);

% u_initial = zeros(np,1);
% for i=1:mnp
%     u_initial(i) = u_0(p(1,i),p(2,i));
% end %for, i

for i=1:np
%     loc2glb = t(1:3,i);
%     x=p(1,loc2glb);
%     y=p(2,loc2glb);
    %xi_old(i) = 1/3 * (u_0(x(1),y(1))+u_0(x(2),y(2))+u_0(x(3),y(3)));
    xi_old(i) = u_0(p(1,i),p(2,i));
end %for, i

size(C)
size(M)
size(xi_old)
size(p)
Xi=xi_old;
for i=1:N
    xi_new = (+C/2+M/delta_t)\((M/delta_t-C/2)*xi_old);


    pdesurf(p,t,xi_new)
    drawnow;
end %for, i


%load vector
b = Assembleb(p,e,t,f);
%replace the rows corresponding to the boundary nodes by corresponding rows
%of I
A(e(1,:),:) = I(e(1,:),:);
%put the boundary value into the RHS
b(e(1,:))=0;

%plot(linspace(1,380,380),A\b)

function C = ConvectionAssembler2D(p,t,bx,by)
    np=size(p,2);
    nt=size(t,2);
    C=sparse(np,np);

    for i=1:nt
        loc2glb = t(1:3,i);
        x=p(1,loc2glb);
        y=p(2,loc2glb);
        [area,b,c] = HatGradients(x,y);
        bxmid = mean(bx(loc2glb));
        bymid = mean(by(loc2glb));

        CK=ones(3,1)*(bxmid*b+bymid*c)'*area/3;
        C(loc2glb,loc2glb) = C(loc2glb,loc2glb) + CK;
    end %for
end %function




function A = AssembleA(p,e,t);
%number of nodes
np = size(p,2);
%number of triangles
nt = size(t,2);
%initialisation
A = sparse(np,np);

%iterating through each triangle
for K = 1:nt
    %actual nodes
    nodes = t(1:3,K);
    %coordinates x and y
    x = p(1,nodes);
    y = p(2,nodes);

    %area
    K_area = polyarea(x,y);

    %gradient components
    bi = (1/2)*[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/K_area;
    ci = (1/2)*[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/K_area;

    %element for the stiffness matrix A
    AK = (bi*bi' + ci*ci') * K_area;

    %add element to A
    A(nodes, nodes) = A(nodes,nodes) + AK;

end %for, K
end %function

function M = MassAssembler2D(p,t)
np = size(p,2); %#nnodes
nt = size(t,2); %#elements
M = sparse(np,np); %allocate mass matrix
for K = 1:nt %loop over elemnts
    loc2glb = t(1:3,K); %local-to-global map
    x = p(1,loc2glb); % node x-coordinates
    y = p(2,loc2glb); % node y
    area = polyarea(x,y); %triangle area
    MK = [2,1,1;1,2,1;1,1,2]/12*area; %element mass matrix
    M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK; %add element masses to M
end %for
end %function



function b = Assembleb(p,e,t,f)
%number of nodes
np = size(p,2);
%number of triangles
nt = size(t,2);
%initialisation
b = sparse(np,1);



%iterating through each triangle
for K = 1:nt
    %actual nodes
    nodes = t(1:3,K);
    %coordinates x and y
    x = p(1,nodes);
    y = p(2,nodes);

    %area
    K_area = polyarea(x,y);
    bK = [f(x(1),y(1));f(x(2),y(2));f(x(3),y(3))]/3*K_area;

    b(nodes) = b(nodes) + bK;
end %for, K
end %function

function [area,b,c] = HatGradients(x,y)
    area = polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end% function

function value = divergence(f, u, x1, x2)
        v = 2*pi*[-x2,x1];
        w = [(u(x1+0.01,x2)-u(x1-0.01,x2))/0.02;(u(x1, x2+0.01)-u(x1,x2-0.01))/0.02]
end %function


% %creating load vector
% function b = LoadAssembler2D(p,t,f)
% np = size(p,2); %#nnodes
% nt = size(t,2); %#elements
% b = zeros(np,1);
%
% for K = 1:nt
%     loc2glb = t(1:3,K); % local-to-global map
%     x = p(1,loc2glb); % node x-coordinates
%     y = p(2,loc2glb); % node y
%     area = polyarea(x,y);
%     bK = [f(x(1),y(1));f(x(2),y(2));f(x(3),y(3))]/3*area; %element load vector
%     b(loc2glb) = b(loc2glb) + bK; %add element loads to b
% end %for
% end %function
%


%
% %continous piecewise linear interpolation
% x = p(1,:); y = p(2,:); %node coordinates
% pif = x.*y; %nodal values of interpolant !!!
%
% %Assembly of the mass matrix M
% M = MassAssembler2D(p,t);
%
% %Assembly of the load vector b
% b = LoadAssembler2D(p,t,f);
%
% %Assembly of the stiffness matrix A
% A = StiffnessAssembler2D(p,t,0);  %??
%
% %Backward Euler
% xi  =   u_0(0.3,0);
%
% for l=1:max(size(time))
%     b_old = zeros(380,1);
%     xi_old = xi(1,:)';
%     xi_new = (M+k*A)\(M*xi_old+k*b_old);
%
%
% end %for, l
% xi_new
%
%
% % %L2Projector2D(f, hmax)
% %
% % nt = size(t,2);
% % Cvel = 0.25;
% % Crv = 1.0;
% % Res = L2Projector2D(f, hmax);
%
%
%
%
% function [A] = Assemble_A(p,t)
% np = size(p,2); nt = size(t,2);
% A = sparse(np,np);
% for K = 1:nt
%     loc2glb = t(1:3,K); % local to global mapping
%     x = p(1,loc2glb); y = p(2,loc2glb);
%     K_area = polyarea(x,y);
%     % compute hat gradient components for K:
%     bi = (1/2)*[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/K_area;
%     ci = (1/2)*[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/K_area;
%     % compute local element stiffness matrix A^k:
%     AK = (bi*bi' + ci*ci') * K_area;
%     % compute global stiffness matrix M:
%     A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK;
% end %for
% end %function
%
%
%
%
%
% %creating Mass Matrix
% function M = MassAssembler2D(p,t)
% np = size(p,2); %#nnodes
% nt = size(t,2); %#elements
% M = sparse(np,np); %allocate mass matrix
% for K = 1:nt %loop over elemnts
%     loc2glb = t(1:3,K); %local-to-global map
%     x = p(1,loc2glb); % node x-coordinates
%     y = p(2,loc2glb); % node y
%     area = polyarea(x,y); %triangle area
%     MK = [2,1,1;1,2,1;1,1,2]/12*area; %element mass matrix
%     M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK; %add element masses to M
% end %for
% end %function
%
%
% function [area,b,c] = HatGradients(x,y)
% area = polyarea(x,y);
% b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
% c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
% end% function
%
% %creating stiffness matrix
% %A = StiffnessAssembler2D(p,e,t)
%
% function A = StiffnessAssembler2D(p,t,a)
% np = size(p,2);
% nt = size(t,2);
% A = sparse(np,np); % allocate stiffness matrix
% for K = 1:nt
%     loc2glb = t(1:3,K); % local-to-global map
%     x = p(1,loc2glb); % node x-coordinates
%     y = p(2,loc2glb); % node y
%     [area,b,c] = HatGradients(x,y);
%     xc = mean(x);
%     yc = mean(y); % element centroid
%     abar = 1;%a(xc,yc); % value of a(x,y) at centroid
%     AK = abar*(b*b'...
%         +c*c')*area; % element stiffness matrix
%     A(loc2glb,loc2glb) = A(loc2glb,loc2glb) ...
%         + AK; % add element stiffnesses to A
% end %for
% end %function
%
%
% %creating load vector
% function b = LoadAssembler2D(p,t,f)
% np = size(p,2); %#nnodes
% nt = size(t,2); %#elements
% b = zeros(np,1);
%
% for K = 1:nt
%     loc2glb = t(1:3,K); % local-to-global map
%     x = p(1,loc2glb); % node x-coordinates
%     y = p(2,loc2glb); % node y
%     area = polyarea(x,y);
%     bK = [f(x(1),y(1));f(x(2),y(2));f(x(3),y(3))]/3*area; %element load vector
%     b(loc2glb) = b(loc2glb) + bK; %add element loads to b
% end %for
% end %function
%
% function Pf = L2Projector2D(f, hmax)
% g = @circleg;
% [p,e,t] = initmesh(g, 'hmax', hmax); %create mesh
% M = MassAssembler2D(p,t); %assemble mass matrix
% b = LoadAssembler2D(p, t, f); %assemble load vector
% Pf = M\b; %solve linear system
% pdesurf(p,t,Pf);
% end %function


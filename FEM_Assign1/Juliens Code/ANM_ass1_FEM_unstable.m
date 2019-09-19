function [ U, M ] = ANM_ass1_FEM_unstable( p, e, t )


%initial data
r_0 = 0.25;
x1_0 = 0.3;
x2_0 = 0;
T = 1;
f = @(x, y) 0;
f_der = @(x1, x2) 2*pi*[-x2,x1]  ;
CFL = 0.5;

%initial function
u_0 = @(x1,x2) 1/2*(1-tanh(((x1-x1_0)^2+(x2-x2_0)^2)/(r_0^2)-1));


%time partition
N = 51;
time = linspace(0, 1, N+1);
delta_t = T/N;


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



for i=1:np
    xi_old(i) = u_0(p(1,i),p(2,i));
    %xi_old(i) = u_0_2(p(1,i),p(2,i), x1_0, x2_0, r_0);
end %for, i

U = zeros(np, N+1);
U(:,1) = xi_old;

for i=1:N

    A1 = C/2+M/delta_t;
    A2 = M/delta_t-C/2;
    bb = A2*xi_old;
    xi_new = A1\bb;
    xi_old = xi_new;
    U(:,i+1)=xi_old;
end %for, i

e = (U(:,end)-U(:,1));


function C = ConvectionAssembler2D(p,t,bx,by)
np=size(p,2);
nt=size(t,2);
C=sparse(np,np);

for i=1:nt
    loc2glb = t(1:3,i);
    x=p(1,loc2glb);
    y=p(2,loc2glb);
    [area,b,c] = HatGradients(x,y);

    bxmid = mean(bx(y));
    bymid = mean(by(x));

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

function value = u_0_2(x1,x2, x1_0, x2_0, r_0)
    if((x1-x1_0)^2+(x2-x2_0)^2<=r_0^2)
        value = 1;
    else
        value = 0;
    end%if
end %function
end


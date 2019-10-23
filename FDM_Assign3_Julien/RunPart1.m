%diagonal matrix
lambda = [-1,0;0,1];

m=129; %problemstorlek
h=1/(m-1);

% Second order accurate SBP operators. Specify nummber of gripdpoints and h
% Note that h depend on the size of the computational domain (it is not
% always 1

%D2=HI(-A+BD)

%% SBP-SAT-approximation
e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

H=(eye(m,m));H(1,1)=0.5;H(m,m)=0.5;
H=h*H;
HI=inv(H);

D1=((.5*diag(ones(m-1,1),1)-.5*diag(ones(m-1,1),-1)));
D1(1,1)=-1;D1(1,2)=1;D1(m,m-1)=-1;D1(m,m)=1;
D1(m,m-1)=-1;D1(m,m)=1;
D1=D1/h;

Q=H*D1 + 1/2*e_1*e_1' - 1/2*e_m*e_m';

D2=((diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)-2*diag(ones(m,1),0)));
D2(1,1)=1;D2(1,2)=-2;D2(1,3)=1;
D2(m,m-2)=1;D2(m,m-1)=-2;D2(m,m)=1;
D2=D2/h^2;

S_U=[-3/2, 2, -1/2]/h;
S_1=zeros(1,m);
S_1(1:3)=S_U;
S_m=zeros(1,m);
S_m(m-2:m)=fliplr(-S_U);

M=-H*D2-e_1*S_1+e_m*S_m;

%% Initial data

eps = 1;
mu = 1;
rr = 0.1; %that's r*

%boundaries
x_l = -1;
x_r = 1;
L = x_r - x_l;
space = linspace(x_l, x_r, m);

%time partition
T = 1.8; %that's t
dt = 0.1 * h;

%Gaussian profiles
theta1 = @(x,t) exp(-((x-t)/rr)^2);
theta2 = @(x,t) -exp(-((x+t)/rr)^2);

%Initial data
u_1 = @(x) theta2(x,0) - theta1(x,0);
u_2 = @(x) theta2(x,0) + theta1(x,0);

%% Solving the ODE
ODEmatrix = kron(lambda,D1);
f = @(t,v) ODEmatrix * v;
v1 = zeros(m,1);
v2 = zeros(m,1);
for i=1:m
    v1(i) = u_1(space(i));
    v2(i) = u_2(space(i));
end %for, i
y_0 = [v1;v2];
y = RungeKutta4(f, 0, y_0, h, m);

for i=1:m
    v1 = y(1:m,i);
    v2 = y(m+1:end,i);
    plot(space, v1, 'r', space, v2, 'b');
    pause(.1)
    drawnow;
end %for, i


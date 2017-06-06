% This code solves the Advection-Diffusion Equation, 
%       u_t+au_x=bu_xx,  
% using a Fraction Step Method (discussed in Chpt 11 of FDM by Leveque).
% This method solves  u_t+au_x=0  over a time step k/2 to get a grid value
% U*.  U* is used in u_t=bu_xx over a second time step k/2 to get the grid
% value, U at time t_0+k.

% Upwind_Advection.m from ``Numerical Schemes'' will be borrowed for the
% advection part of this solver, and Two_Point_BVP_Centered_Difference.m
% for the diffusion.

a=0.0475;  % advection coefficient
b=1e-2;    % diffusion coefficient

f=@(x,t) sin(8*pi*x);

m=100;   % # spatial grid points
xend=1;    % final x value
h=xend/m;   % step size
x=(0:h:xend);   % x grid

n=1000;
tend=2;
k=tend/n;
t=(0:k:tend);

if (abs(a*k/h)>1)
	error(['Please select "a" smaller than ', num2str(h/k), '. ',...
        num2str(0.95*h/k), ' gives |ak/h|=0.95.']) % Upwind stable for |ak<h|<1
end

% Setup Centered Difference method for Diffusion

A=sparse(m+1,m+1);   % Sparse uses less memory than zeros
for i=2:m;       % Rows 1 and m+1 reserved for BC
    A(i,i)=2;    % main diagonal is 2
    A(i,i-1)=-1; % upper diagonal is -1
    A(i,i+1)=-1; % lower diagonal is -1
end
A=h^(-2)*A;

C=sparse(eye(m+1,m+1)+k*b*A);   % If F(1)=0, Dirichlet BC satisfied
C(m+1,m)=-1/h; C(m+1,m+1)=1/h;  % Neumann BC at x=1 (we will set F(m+1)=0)

% Advection, then Diffusion

u1=zeros(1,length(x));
u1(round(0.1*(length(x)/xend)+1):round(0.3*(length(x)/xend)))=1;
IC=u1;  % Save initial condition for plotting purposes

close all
figure

for i=1:n
	u0=u1;	% update (i-1)k solution vector
    if a>0
	u1(2:m)=u0(2:m)-(a*k/h)*(u0(2:m)-u0(1:m-1));  % upwind for a>0
    else
	u1(2:m)=u0(2:m)-(a*k/h)*(u0(3:m+1)-u0(2:m));  % upwind for a<0
    end
	u1(1)=0; 	% Dirichlet BC at x=0
    u1(m+1)=u1(m);  % Neumann BC at x=1
    F=zeros(m+1,1);   
    for j=2:m
        F(j)=k*f(x(j),t(i))+u0(j);   
    end
    u1=C\F; % centered difference for diffusion
    plot(x,u1); axis([0 1 -1 1]);
    pause(0.01)
end
plot(x,u1,x,IC)


% This code solves the Advection-Diffusion Equation, 
%       u_t+au_x=bu_xx
% using a Fraction Step Method (discussed in Chpt 11 of FDM by Leveque).
% This method solves  u_t+au_x=0  over a time step k/2 to get a grid value
% U*.  U* is used in u_t=bu_xx over a second time step k/2 to get the grid
% value, U at time t_0+k.

% Upwind_Advection.m from ``Numerical Schemes'' will be borrowed for the
% advection part of this solver, and Two_Point_BVP_Centered_Difference.m
% for the diffusion.

a=0.2;  % advection coefficient
b=1;  % diffusion coefficient

I=@(y) -(y-0.7808).^2+0.6096;  % u(x,0)

m=40;   % # spatial grid points
xend=1;    % final x value
h=xend/m;   % step size
x=[0:h:xend];   % x grid

n=20;
tend=2;
k=tend/n;
t=[0:k:tend];

if (abs(a*k/h)>1)
	error(['Please select "a" smaller than ', num2str(h/k), '. ',...
        num2str(0.95*h/k), ' gives |ak/h|=0.95.']) % Upwind stable for |ak<h|<1
end

% Setup Centered Difference method for Diffusion

A=zeros(m,m);   % Matrix for numerical scheme.
A(1,1)=2;  % Columns 1 and m done manually because they do not
A(1,2)=-1;   %   include both the lower and upper diagonals.
A(m,m)=h;  
A(m,m-1)=-h;
for i=2:m-1; 
    A(i,i)=2;  % main diagonal is -2
    A(i,i-1)=-1; % upper diagonal is 1
    A(i,i+1)=-1; % lower diagonal is 1
end
invA=inv(A);

% Advection, then Diffusion

u1=I(x);  % initial data u(x,0)

close all
figure
plot(x,u1)
hold on
for i=1:n
	u0=u1;	% update (i-1)k solution vector
	u1(2:m)=u0(2:m)-(a*k/h)*(u0(2:m)-u0(1:m-1)); % upwind advection for a>0
	u1(1)=0; 	% Dirichlet BC at x=0
    u1(m+1)=u1(m);  % Neumann BC at x=1
    F=zeros(1,m);
    for i=2:m
        F(i-1)=(u1(i)-u0(i))/k;   % f(x(i))=u_t(x(i))
    end
    F(m)=0;
    u1(2:m+1)=(2.1*h^2.*invA*F')'; % centered difference for diffusion
    plot(x,u1)
    pause(0.1)
end
hold off

% Multiple species advection diffusion reaction (3 species)

% This code solves the Advection-Diffusion Equation, 
%       u_t+au_x=bu_xx,  
% using a Fraction Step Method (discussed in Chpt 11 of FDM by Leveque).
% This method solves  u_t+au_x=0  over a time step k/2 to get a grid value
% U*.  U* is used in u_t=bu_xx over a second time step k/2 to get the grid
% value, U at time t_0+k.

% Upwind_Advection.m from ``Numerical Schemes'' will be borrowed for the
% advection part of this solver, and Two_Point_BVP_Centered_Difference.m
% for the diffusion.

a=0.00475;  % advection coefficient
b=1e-2;    % diffusion coefficient
K1=0.1;  % reaction coefficientn
K2=0.001;   % reaction coefficient
dBC=2;  % Dirichlet BC at x=0
nBC=-1;  % Nuemann BC at x=1 (or xend, if different than 1)

tol=1e-8; % tolerance for exiting Multi_Newton loop
iter=20; % number of iterations for Multi_Newton

m=100;   % # spatial grid points
xend=1;    % final x value
h=xend/m;   % step size
x=(0:h:xend);   % x grid

n=1000;
tend=1;
k=tend/n;
t=(0:k:tend);

figure

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

% Initial Profiles 

u1=ones(1,length(x));
u1(round(0.1*(length(x)/xend)+1):round(0.3*(length(x)/xend)))=1;

u2=ones(1,length(x));
u2(round(0.1*(length(x)/xend)+1):round(0.4*(length(x)/xend)))=1;

u3=ones(1,length(x));
u3(round(0.1*(length(x)/xend)+1):round(0.5*(length(x)/xend)))=1;

u=[u1; u2; u3];

u=rand(3,length(x));
uL=u;

fun1=@(x) (x-0.5).^2+dBC-0.25;
fun2=@(x) x+dBC;
fun3=@(x) x.^4+dBC;

u=[fun1(x);fun2(x);fun3(x)];
uL=u;

figure

f=@(x,t) 0;% sin(8*pi*x);

time=t(1);
i=0;    % keep track of number of iterations in while loop

for i=1:n
%while time<tend
    
    i=i+1;
    
    %_________________________________________________________
    
    % Advection and Reaction first
    
    u=advection(a,k,dBC,nBC,h,m,u,3);  % BC prescribed in "advection.m"
    [u,num]=Multi_Newton(k,m,u,tol,iter,K1,K2);
    
    %__________________________________________________________
    
    % Complete diffusion step
    
    for j=1:3
    F=zeros(m+1,1);
    F(1)=dBC;  % Dirichlet BC
    F(m+1)=nBC; % Neumann BC
    for l=2:m
        F(l)=k*f(x(l),t(i))+u(j,l);   
    end
    u(j,:)=C\F;     % centered difference for diffusion
    end
    
    %___________________________________________________________
    
    % Update time step based on number of iterations it takes for Newton to
    % converge within given tolerance; advance in time accordingly
    %%{
    if num<4
        k=k*2;
    else if num>8
            k=k/2;
        end
    end
    %}
    time=time+k;
    K(i)=k;     % Keep track of changes to time step

    %___________________________________________________________
    
    plot(x,u(1,:),'r',x,u(2,:),'b',x,u(3,:),'k',x,uL(1,:),'r--',x,uL(2,:),'b--',x,uL(3,:),'k--')
    %plot(x,u(1,:),'r',x,u(2,:),'b',x,u(3,:),'k');
    pause(0.01)
end
K

plot(x,u(1,:),'r',x,u(2,:),'b',x,u(3,:),'k',x,uL(1,:),'r--',x,uL(2,:),'b--',x,uL(3,:),'k--')
xlabel('Spatial Position')
ylabel('Concentration')
lgn=legend('A solution','B Solution','C Solution','A Initial','B Initial','C Initial');
lgn.Location='northwest';
function out=Multi_Newton_Test_1funct(k,u,tol,N)

clear uG

% k := time mesh size
% m := # of spatial grid points
% u := concentration after advective step
% tol := accuracy requirement for exiting Newton loop
% N := max number of Newton iterations

close all
figure
hold on

out=NaN;

% Reaction equations

K1=1;
K2=1;

f1=@(x) -K1*x(1,1)*x(2,1)+K2*x(3,1);
f2=@(x) -K1*x(1,1)*x(2,1)+K2*x(3,1);
f3=@(x) K1*x(1,1)*x(2,1)-K2*x(3,1);

% Guesses for u(1,:),u(2,:),u(3,:)

uG(1,1)=u(1,1)+0.01;
uG(2,1)=u(2,1)+0.01;
uG(3,1)=u(3,1)+0.01;

% Initalize residual

res=zeros(size(u));

% Solve for updated u value

for i=1:N
    
    res(1,1)=uG(1,1)-u(1,1)-k.*f1(uG);
    res(2,1)=uG(2,1)-u(2,1)-k.*f2(uG);
    res(3,1)=uG(3,1)-u(3,1)-k.*f3(uG);
    
    if norm(res,inf)<tol
        out=uG;
        i
        return
    end
    
    %Jac=zeros(3,3); % Jacobian matrix

    Jac11=1-k*(-K1*uG(2,1));     % dr1/duG1
    Jac21=-k*(-K1*uG(2,1));      % dr2/duG1
    Jac31=-k*(K1*uG(2,1));       % dr3/duG1
    
    Jac12=-k*(-K1*uG(1,1));      % dr1/duG2
    Jac22=1-k*(-K1*uG(1,1));     % dr2/duG2
    Jac32=-k*K1*uG(1,1);         % dr3/duG2
    
    Jac13=-k*K2;             % dr1/duG3
    Jac23=-k*K2;             % dr2/duG3
    Jac33=1-k*(-K2);         % dr3/duG3
    
    % Update guess
    
    Jac=[Jac11,Jac12,Jac13;Jac21,Jac22,Jac23;Jac31,Jac32,Jac33];
    step=-Jac\res;
    uG(:,1)=uG(:,1)+step;
    
    plot([1,2,3],uG(:,1)')
    pause(0.1)
end

hold off

end

function [out,num]=Multi_Newton(k,m,u,tol,N,K1,K2)

%clear uG

out=NaN;
num=NaN;

% k := time mesh size
% m := # of spatial grid points
% u := concentration after advective step
% tol := accuracy requirement for exiting Newton loop
% N := max number of Newton iterations
% K1, K2 := reaction rate constants

% Reaction equations

f1=@(x) -K1.*x(1,:)*(x(2,:)')+K2.*x(3,:);
f2=@(x) -K1.*x(1,:)*(x(2,:)')+K2.*x(3,:);
f3=@(x) K1.*x(1,:)*(x(2,:)')-K2.*x(3,:);

% Guesses for u(1,:),u(2,:),u(3,:)

uG(1,:)=u(1,:)+0.01*ones(size(u(1,:)));
uG(2,:)=u(2,:)+0.01*ones(size(u(2,:)));
uG(3,:)=u(3,:)+0.01*ones(size(u(3,:)));

% Initalize residual

res=zeros(size(u));

% Solve for updated u value

for i=1:N
    
    res(1,:)=uG(1,:)-u(1,:)-k.*f1(uG);
    res(2,:)=uG(2,:)-u(2,:)-k.*f2(uG);
    res(3,:)=uG(3,:)-u(3,:)-k.*f3(uG);
    
    
    % Tolerance applied to Actual Error
    %%{
    if norm(res,inf)<tol
        out=uG;
        num=i-1;
        return
    end
    %%}
    %{
    % Tolerance applied to Relative Error
    
    Relative_res=res./uG;
    
    if norm(Relative_res)<tol
        out=uG;
        num=i-1;
        return
    end
    %}

    Jac11=ones(size(uG(1,:)))-k.*(-K1.*uG(2,:));     % dr1/duG1
    Jac21=-k.*(-K1.*uG(2,:));                        % dr2/duG1
    Jac31=-k.*(K1.*uG(2,:));                         % dr3/duG1
    
    Jac12=-k.*(-K1.*uG(1,:));                        % dr1/duG2
    Jac22=ones(size(uG(2,:)))-k.*(-K1.*uG(1,:));     % dr2/duG2
    Jac32=-k.*(K1.*uG(1,:));                         % dr3/duG2
    
    Jac13=-k.*(K2.*ones(size(uG(3,:))));             % dr1/duG3
    Jac23=-k.*(K2.*ones(size(uG(3,:))));             % dr2/duG3
    Jac33=ones(size(uG(3,:)))-k.*(-K2.*ones(size(uG(3,:))));  % dr3/duG3
    
    % Update guess
    
    for j=1:m+1
        Jac=[Jac11(j),Jac12(j),Jac13(j);Jac21(j),Jac22(j),Jac23(j);Jac31(j),Jac32(j),Jac33(j)];
        step=-Jac\res(:,j);
        uG(:,j)=uG(:,j)+step;
    end
    
    %plot([1,2,3],uG(:,1)','r',[1,2,3],uG(:,2)','b')
    %pause(0.05)
end

end

%}
  
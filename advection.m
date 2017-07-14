% Advective Step

function u=advection(a,k,dBC,nBC,h,m,u0,N)
 
 % a := advection coefficient
 % k := time mesh size
 % dBC := Dirichlet BC (x=0)
 % nBC := Neumann BC (x=1)
 % h := spacial mesh size
 % m := # of spatial grid points
 % u0 := chemical profile to be advected
 % N := number of chemical species, u(1,:),...,u(N,:)

for j=1:N;
   if a>0
    u(j,2:m+1)=u0(j,2:m+1)-(a*k/h)*(u0(j,2:m+1)-u0(j,1:m));  % upwind for a>0
    u(j,1)=dBC; % advection uses 1 BC at the inflow boundary
   else
    u(j,1:m)=u0(j,1:m)-(a*k/h)*(u0(j,2:m+1)-u0(j,1:m));  % upwind for a<0
    u(j,m+1)=u(j,m)+h*nBC;  % a<0 => inflow boundary at x=1
   end
end
 
end

 
 
 
 
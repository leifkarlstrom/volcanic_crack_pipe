function [rho K c a b, p, n]=magma_st(x,BG)
%compute the depth dependent magma static properties given depth x 
%using ode45. All the parameters related to the equation 
%of state of magma is built into magma_eos.m.
%Input:
%       x:  the grid points (Note that x=0 means surface and increases as depth increases)
%
%Output:
%       n = gas mass fraction
%       rho = density
%       c = sound speed
%       K = bulk modulus
%       a = gas exsolution, -(drho/dn)/rho
%       b = gas exsolution, -dneq/dp
%
%Note: to be constant with the sign notation of the rest of the code.
%first entry means the bottom while the last entry means the surface.
%Refer magma_eos.m for more parameters related to the eos of magma.

% Integrate dp/dz=rho*g with rho=rho(p,n_eq);

% %gravitational acceleration.
% g=10;
 %

%reference pressure-atmosphere pressure in this case
 p_ref=1e5; %Pa
% 
% %number of grid points
% N=length(x);
% 
% %allocate storage for initial pressure
% p=zeros(size(x));
% 
% % Numeric integration using Forward Euler.
% 
% p(N)=p_ref;
% for i=N:-1:2
%     
%     rho_i=magma_eos(p(i));
%     
%     dx_i=x(i)-x(i-1);
%     
%     Dp=g*rho_i*dx_i;
%     
%     p(i-1)=p(i)+Dp;
% end

%Numeric integration using ode45.
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
sol = ode45(@(x,p) f(x,p,BG),[x(end) x(1)],p_ref,options);
%evaluate solution at high accuracy
p=deval(sol,x)';

    if strcmp(BG,'specified_n')
        [rho,~,~, K, c, a, b,n]=magma_eos(p,BG,x);
    else
        [rho,~,~, K, c, a, b,n]=magma_eos(p,BG);
    end

%[rho,~,~, K, c, a, b,n]=magma_eos(p,BG);
%keyboard
end

% define function dp/dx.
    function f=f(x,p,BG)
        g=9.8; %gravity m/s2

        if strcmp(BG,'specified_n')
            rho=magma_eos(p,BG,x);
        else
            rho=magma_eos(p,BG);
        end

        f=-rho*g;
    end


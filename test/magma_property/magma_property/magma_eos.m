function [rho, rho_g, rho_l, K, c, a, b,n_eq]=magma_eos(p,BG,z,n)
% magma_eos(p,n) returns the general values whil magma_eos(p) return
% equibrium values. 
% The optional input z is for prescribing gas content in one BG option
%
% p: pressure
% n: gas mass fraction
% rho = density of mixture
% rho_g: gas density
% rho_l: liquid density
% c = sound speed
% K = bulk modulus
% a = gas exsolution, -(drho/dn)/rho
% b = gas exsolution, -dneq/dp

% equation of state:
% a mixture of ideal gas and melt with constant compressibility.
% gas: rho_g=p/(R*T_m).
% melt:rho_l=rho_l_ref*(1+(p-p)*beta_l).

n_t =3e-3;% total volative content 0.05% w.pct of water

R = 461.5; %[J/kg/K] specific gas constant for water vapor
T_m =1155+273.15;% Temperature of magma about 1000 deg C.

beta_l = 1e-10; %[Pa^-1] compressibility constant for magma melt
K_l = 1/beta_l;% [Pa] bulk modulus of magma melt
rho_l_ref = 2900;%[kg/m^3], reference density, typical density for basalt
p_ref = 1e5;%[Pa], reference pressure, atmosphere pressure

% Constants for Henry's law. n_eq=n_o-s*p^m;
s_h = 4e-6;%[Pa^(-1/2)]. solubility constant.
m = 0.5;% Henry's law is 1/2 for water, would be 1 for CO2

switch BG
    
    case 'henrys_law'
p_ex = (n_t/s_h)^(1/m);%[Pa] exsolution pressure.        
if nargin<3 % if magma_eos(p) return equilibrium value
    n=n_t-s_h*p.^m; %Henry's law
    n(p>p_ex)=0;
    n_eq=n;
end        
%compute density
rho_g = p./(R*T_m);

    case 'Halemaumau_max'        
        %joint equilibrium solubilty of CO2 and H2O from Papale 2006,
        %assuming 0.36 wt% water and 500 ppm CO2 (max esimates from Edmonds
        %et al 2013 for Halemaumau - H2O - and kilauea generally for CO2)
        load('KilaueaMax.mat') 
        load('RK_H2O_CO2_density.mat') 
p_ex = gasfrac(end,1)*1e5;%[Pa] exsolution pressure.        
if nargin<3 % if magma_eos(p) return equilibrium value
    if p>p_ex
    n(p>p_ex)=0;
    else
    n(p<=p_ex) = interp1(gasfrac(:,1)*1e5,gasfrac(:,2),p(p<=p_ex))/100; 
    end
    n=n';
    n_eq=n;
end          

rho_g = interp1(Pgas*1e5,rhogRK,p);

%p./(R*T_m);        
    case 'Halemaumau_mean'
        %joint equilibrium solubilty of CO2 and H2O from Papale 2006,
        %assuming 0.36 wt% water and 340 ppm CO2 (mean esimates from Edmonds
        %et al 2013 for Halemaumau - melt inclusions)
        load('KilaueaMean.mat')
        gasfrac(:,1)=gasfrac(:,1)/1e5;
        load('RK_H2O_CO2_density.mat')
p_ex = gasfrac(end,1)*1e5;%[Pa] exsolution pressure.                
if nargin<3 % if magma_eos(p) return equilibrium value
    if p>p_ex
    n(p>p_ex)=0;
    else
    n(p<=p_ex) = interp1(gasfrac(:,1)*1e5,gasfrac(:,2),p(p<=p_ex))/100; 
    end
    n=n';
    n_eq=n;
end  
%rho_g = p./(R*T_m); 
rho_g = interp1(Pgas*1e5,rhogRK,p);

    case 'specified_n'
        %prescribe the spatial variation in gas mass fraction based on
        %two segment (four points total)
if nargin<4 % if magma_eos(p) return equilibrium value and we have z value
%keyboard
z1 = 700;
z2 = 300;

n1 = 0.004;
n2 = 0.001;

if length(z)==1
    if z>z1
    n=n1;
    elseif z<=z1 & z>z2
    n = n1 + (n2-n1)*(z-z1)/(z2-z1); %linear interpolation 
    else
    n = n2;    
    end
else
    %keyboard
    n = zeros(size(z));
    n(z>z1) = n1;
    n(z<=z1 & z>z2) = n1 + (n2-n1)*(z(z<=z1 & z>z2)-z1)/(z2-z1);
    n(z<=z2) = n2;
end

    n=n';
    n_eq=n;

rho_g = p./(R*T_m); 
end          

end


rho_l = rho_l_ref+(p-p_ref)*beta_l;
rho = 1./(n./rho_g+(1-n)./rho_l);

K = 1./(rho.*(n./(rho_g.*p)+(1-n)./(rho_l*K_l)));
a = rho.*(1./rho_g-1./rho_l);
b = s_h*m*p.^(m-1);
% if p>p_ex, b=0;
%b(p>p_ex)=0;
b(n==0)=0;

%instantaneous wave speed
c=sqrt(K./rho);

end 

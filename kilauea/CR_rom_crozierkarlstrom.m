function [CRout] = CR_rom_crozierkarlstrom(out)
%evaluate natural frequency and quality factor of conduit-reservoir mode
%based on the magmastatic model of Crozier and Karlstrom 2022 SciAdv

Pm.H = out.M.L; %length of the conduit 

% Pm.dz = out.M.dz;
% Pm.rho = out.M.rho;%1000*ones(1,N+1); %density as fctn of depth z kg/m3

N = 1e4; %number of points desired in zvec
Pm.dz = Pm.H/(N -1) ;
z = linspace(0,Pm.H,N); %generate z vec with custom dz

if out.M.interface_split==true 
    %remove the repeated index if using an interface
    %this may be slightly inaccurate if there is a true jump across the
    %interface, but otherwise will be accurate
    out.z(out.M.split_index)=[];
    out.M.mu(out.M.split_index)=[];
    out.p(out.M.split_index)=[];
    out.M.rho(out.M.split_index)=[];
end


oldz = out.z; 

Pm.rho = interp1(oldz,out.M.rho,z);

Pm.drhodz = gradient(Pm.rho,Pm.dz); %density gradient 
Pm.mu = interp1(oldz,out.M.mu,z);%viscosity as fctn of depth z in Pas

%parameter definitions for conduit (need to be consistent with simulation!)
Pm.g = out.M.g; %gravity

Pm.rad = sqrt(out.M.S/pi); %use the conduit cross section to get radius
%out.M.R; %conduit radius as fct of z based on input

%Pm.rad = out.M.R*ones(size(z));%interp1(oldz,out.M.R,z); %interpolate to new zvec

Pm.Ct = out.M.Ct;%total storativity of reservoir (magma+elastic)

Pm.theta = 90; %dip angle (assumed constant in z) of conduit

%evaluate the constants needed to define properties of damped harmonic
%oscillator (the reduced order model of Crozier and Karlstrom 2022)

try
    [c1,c2,c3,Omega] = evaluateparams(Pm);

    %note! angular frequency is related to frequency by factor of 2pi!
    Tvis = 2*pi/Omega;
    
    DeltaRho = Pm.rho(end)-Pm.rho(1); %density at top of column minus base
    %average density 
    rhoavg = 1/Pm.H*trapz(Pm.dz,Pm.rho);%0.5*(Pm.rho(1)+Pm.rho(end));%mean(Pm.rho);
    %inviscid resonant period
    % Tinvis = 2*pi*sqrt(Pm.H*mean(Pm.rho)/ ...
    %     ((Pm.rho(end)-DeltaRho)*Pm.g+pi*Pm.rad(1)^2/Pm.Ct));
    Tinvis = 2*pi*sqrt(Pm.H*rhoavg/ ...
        ((Pm.rho(end)-DeltaRho)*Pm.g+pi*Pm.rad(1)^2/Pm.Ct));
    
    %fully developed resonant period
    Tfd = 2*pi./sqrt(((Pm.rho(end)-DeltaRho)*Pm.g+pi*Pm.rad(1)^2/Pm.Ct)./(Pm.H*rhoavg) - ...
        16*mean(Pm.mu).^2/(Pm.rad(1)^4*rhoavg^2));
    
    lambda = 4*mean(Pm.mu)/(Pm.rad(1)^2*rhoavg);
    
    %quality factor
    Q = Omega*c1/c2;
    
    %fully developed Quality factor
    Qfd = Pm.rad(1)^2 * rhoavg/(8*mean(Pm.mu)) * ...
        sqrt((Pm.g*(Pm.rho(end)-DeltaRho)+ pi*Pm.rad(1)^2/Pm.Ct)/(Pm.H*rhoavg) - ...
        16*mean(Pm.mu)^2/(Pm.rad(1)^4*rhoavg^2));
    
    disp(['Conduit-Reservoir period (ROM) is ' num2str(Tvis) ' sec'])
    disp(['CR Quality factor (ROM) is ' num2str(Q)])
    
    disp(['For reference, inviscid resonant period from Liang is ' num2str(Tinvis) ' sec'])
    disp(['For reference, fully developed flow resonant period is ' num2str(Tfd) ' sec'])
    % disp(['For reference, fully developed Q is ' num2str(Qfd) ])
    
    CRout.T = Tvis;
    CRout.Q = Q;
catch
    CRout.T = "Overdamped";
    CRout.Q = "Overdamped";
    disp('System is overdamped')
end

end

%------------------------functions called----------
function [c1,c2,c3,Omega] = evaluateparams(Pm)

%implents CrozierKarlstrom reduced order model

%inertial term
c1 = Pm.rad(1).^2*trapz(Pm.dz,1./sind(Pm.theta).*Pm.rho./Pm.rad.^2);


%restoring force term
c3 = -Pm.rad(1).^2*(Pm.g*...
    (trapz(Pm.dz,Pm.drhodz./Pm.rad.^2)-(Pm.rho(end)./Pm.rad(end).^2).*sind(Pm.theta)) ...
    - pi/Pm.Ct .*sind(Pm.theta));

%solve for the natural angular frequency of the system
z0=0.03; %starting guess for minimization

%solve implicitly for conduit-reservoir natural angular frequency (inverse of
%period)
if any(Pm.mu==0)
    c2 = NaN;
else
    c2 = dampingterm(z0,Pm);
end

%solve implicitly for resonance including frequency dependent damping
if c2==0||isnan(c2)
   
    Omega = sqrt(c3/c1);
    
else

    Omega = fzero(@(omega) myfun(omega,Pm,c1,c3),z0);
    
end

c2 = dampingterm(Omega,Pm);

end

function c2 = dampingterm(omega,Pm)
%calculate viscous damping term
alpha = sqrt(omega.*Pm.rho./Pm.mu)*1i^(3/2);
integrand = 1./sind(Pm.theta).*(Pm.mu.*alpha)./Pm.rad.^3 ...
    .*besselj(1,Pm.rad.*alpha)./besselj(2,Pm.rad.*alpha);
%keyboard
c2 = 2*Pm.rad(1).^2 * real(trapz(Pm.dz,integrand));
end

function f = myfun(omega,Pm,c1,c3)

c2 = dampingterm(omega,Pm);

%throw an error if overdamped (too high viscosity)
    if c3/c1<(c2/(2*c1))^2%abs(imag(f))>0
        disp(c3/c1)
        disp((c2/(2*c1))^2)
    
        error('System is overdamped, no resonance') 
    else
        %now set up the function to be minimized (move terms in eqn 28 all to one
        %side)
        f = omega - sqrt(c3/c1 - (c2/(2*c1)).^2);

    end


end








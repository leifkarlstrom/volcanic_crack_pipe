%script for running 2d conduit code with quasistatic reservoir at bottom

%based on code Liang et al., 
% 
%added post processing workflow for spectral analysis and plotting and
%comparison to analytic formulas

clear 
close all

%add directories from Liang code to path
source = '../mains';
addpath(genpath(source));
source = '../source';
addpath(genpath(source));

% specify background state 
bgstate = 'parameterized'; %LK note: not using this right now

%call script to especify conduit parameters, based on BGstate
[Mc] = setparameters(bgstate);

%build the model 
Model = conduit_internal_g(Mc);

skip = 1; %only save output every "skip" steps to save memory

T = 250; %total time in sec

eigmodeonly = 1; %look only at eigenmodes, or do full timestepping

% time stepping
use_imex = true; % a flag to choose if to use IMEX or not.
CFL = 0.5;
hmin = min([Model.geom.dz]);
cmax = max(Model.M.c);
dt = CFL*hmin/cmax;
nt = ceil(T/dt); %total number of time steps
time  = [skip:skip:nt]*dt;


if ~ use_imex
    A = Model.Ae + Model.Ai;
    % this matrix A is what you need for analyzing the eigenmode and
    % energetics.
else
    [L,U,p,q,B] = imex_ark4_get_lu(Model.Ai,dt);
    A = Model.Ae;
end

% Storage arrays
%ICs
out.p(:,1) = zeros((Mc.nz+1),1);
out.t(1) = 0; 
out.z = Mc.z;
out.dt = dt;
out.skip = skip;

out.M = Mc;


if eigmodeonly
    %if we are just looking at the resonant T and Q, we can skip the
    %timestepping and just look at eigenvalues of the RHS

[evec,e] = eig(full(A));
e = diag(e);

%find eigenvalues that match target range of imag and real part
mask = abs(imag(e))<20&real(e)>-5&abs(imag(e))>5e-2;
LF = find(mask);

T = 2*pi./imag(e(LF))
Q = abs(imag(e(LF))./(2.*real(e(LF))))

figure
%hold on
plot(real(e),imag(e),'o')
ylabel('imaginary part of (s)')
xlabel('real part of (s)')
%keyboard

else

fun = @(u,t) A*u + Model.Fp(:,1)*Model.M.G(t);
tic

for i=1:nt
    t = (i-1)*dt;
    if ~use_imex
        Model=Model.update(lsrk4(fun,Model.u,t,dt));
    else
        Model=Model.update(imex_ark4_lu(Model.u,t,dt,fun, Model.Ai,L,U,p,q));
    end
    
    if mod(i, round(nt/100))==0
        fprintf( '%% %f  finished', round(i*100/nt));
        toc;
    end
    
    if mod(i,skip) == 0
        iter = i/skip;
        %plot solution.
        nr  = Model.geom.nr;
        nz = Model.geom.nz;
        
        % velocity
        vz  = reshape(Model.field(Model.u, 1), [nr,nz]); % [nr, nz]
        
        % pressure
        pz  = Model.field(Model.u, 2); % [nz, 1]
        
        % displacement
        h    = Model.field(Model.u, 3); % [nz,  1]

        % lake height
        hL    = Model.field(Model.u, 4); % [nz,  1]
        
        % width averaged velocity.
        uz  = (Model.op.W1*vz)';%   [nz, 1]

        if strcmp(Mc.BCtype,'quasistatic')
            p_c = Model.field(Model.u, 5);

            out.p_c(:,iter+1) = p_c;
        end
        
        z   = Model.geom.z;
        rm  = Model.geom.rm;
        
        time_i = i*dt;

        out.p(:,iter+1) = pz;
        out.h(iter+1) = hL;
        out.t = [out.t, time_i];
       
    end
end

%% Now we have run the model, now do post-processing

%if strcmp(Mc.BCtype,'quasistatic')
[Fv,FTs,Iv,spectrum] = compute_fft(out.p_c,out.dt);
%end
%in conduit 
[Fv2,FTs2,Iv2,spectrum2] = compute_fft(out.p(2,:),out.dt);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1
periods = 1./Fv2;
out.periods = periods;
out.spectrum = spectrum;
out.spectrum2 = spectrum2;

plotsolutionfields(out)

% Identify peaks in spectrum
%[peaks, peak_locs]=findpeaks(spectrum,'MinPeakProminence',1.5);
[peaks, peak_locs] = findpeaks(spectrum2, 'MinPeakHeight', 5); % Adjust threshold as needed
disp(['peaks in spectrum ' num2str(periods(peak_locs)) ' sec'])

end

[CRout] = CR_rom_crozierkarlstrom(out);

%this implements eqn S36 in Crozier/Karlstrom 
% CRmode_constR = 2*pi*sqrt(out.M.L*0.5*(out.M.rho(end)+out.M.rho(1))/...
%     (out.M.g*(out.M.rho(end)-(out.M.rho(end)-out.M.rho(1))) + pi*out.M.R(1)^2/out.Ct));

%but version with mean of full density seems more accurate
% CRmode_constR = 2*pi*sqrt(out.M.L*mean(out.M.rho)/...
%     (out.M.g*(out.M.rho(end)-(out.M.rho(end)-out.M.rho(1))) + pi*out.M.R(1)^2/out.Ct));

OrganPipe_OO = 2*out.M.L/mean(out.M.c);

%CRmode_constR
disp(['Organ Pipe open-open 1st mode based on mean c is ' num2str(OrganPipe_OO) ' sec'])


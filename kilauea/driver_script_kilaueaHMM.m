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
source = 'wilde_magmastatic';
addpath(genpath(source));
source = 'MELTS_lookup_tables';
addpath(genpath(source));



%call script to specify conduit parameters, based on BGstate
[Mc] = setparameters();
keyboard



%build the model 

% if interface is split, conduit radius will set the resolution for full
% column -KW
% Mc.R = Mc.Rcond;
Model = conduit_internal_g(Mc);

skip = 1; %only save output every "skip" steps to save memory

T = 200; %total time in sec

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
out.p(:,1) = zeros(Model.geom.nz,1);%zeros((Mc.nz+1),1);

out.t(1) = 0; 
out.z = Model.geom.z;%Mc.z;
out.dt = dt;
out.skip = skip;

out.M = Mc;

if eigmodeonly
    %if we are just looking at the resonant T and Q, we can skip the
    %timestepping and just look at eigenvalues of the RHS

% evec has all of the spatial information (these are the eigenvectors)
[evec,e] = eig(full(Model.Ae + Model.Ai));
e = diag(e);
if isempty(e)
    warning('No eignevalues.');
end

%find eigenvalues that match target range of imag and real part
mask = abs(2*pi./imag(e))>.5 & real(e)>-.2 & abs(imag(e))>5e-2;
%abs(imag(e))<20 & real(e)>-5 & abs(imag(e))>5e-2;
LF = find(mask);
% evec = evec( : , LF);

if isempty(imag(e(LF)))
    warning('No eignevalues in range.');
end

T = 2*pi./imag(e(LF));
Q = abs(imag(e(LF))./(2.*real(e(LF))));

% just keep the positive ones
% mask = T>0;
% indices = find(mask);
% 
% T = T(indices);
% Q = Q(indices);
% evec = evec(: , indices);
display(T);
display(Q);

figure
%hold on
plot(real(e),imag(e),'o')
ylabel('imaginary part of (s)')
xlabel('real part of (s)')
%keyboard

%extract dimensions of the model variables: 
% [velocity (not xsectionally avg), pressure, displacement, surface height, chamber pressure]
% [vz, pz, h, hL, p_c]
Dims = Model.dimensions(); 

for i=1:length(T)
    lbl{i} = ['T = ' num2str(round(T(i)*10)/10) ', Q = ' num2str(round(Q(i)*10)/10)]; 
end

% This plot is giving an error, need to fix #KW
% LF exceeds the array limit
figure
subplot(2,1,1)
%plot(Model.M.z,real(evec(Dims(1)+1:Dims(1)+Dims(2),LF)));
plot(Model.geom.z,real(evec(Dims(1)+1:Dims(1)+Dims(2),LF)));
ylabel('Re[pressure eig vec]')
xlabel('distance (m)')
subplot(2,1,2)
%plot(Model.M.z,imag(evec(Dims(1)+1:Dims(1)+Dims(2),LF)));
plot(Model.geom.z,imag(evec(Dims(1)+1:Dims(1)+Dims(2),LF)));
ylabel('Im[pressure eig vec]')
xlabel('distance (m)')
legend(lbl)

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
if strcmp(Mc.BCtype,'quasistatic')
[Fv,FTs,Iv,spectrum] = compute_fft(out.p_c,out.dt);
end
%in conduit 
[Fv2,FTs2,Iv2,spectrum2] = compute_fft(out.p(2,:),out.dt);

%Define the frequency domain f and plot the single-sided amplitude spectrum P1
periods = 1./Fv;
out.periods = periods;
out.spectrum = spectrum;
out.spectrum2 = spectrum2;

plotsolutionfields(out)


% Make T and Q estimations based on model output
[T_condres, Q_condres, T_ac, Q_ac] = post_process(out.t,out.p_c,out.dt);

% Identify peaks in spectrum
% [peaks, peak_locs]=findpeaks(spectrum,'MinPeakProminence',1.5);
% [peaks, peak_locs] = findpeaks(spectrum2, 'MinPeakHeight', 5); % Adjust threshold as needed
% [peaks, peak_locs] = findpeaks(spectrum, 'SortStr','descend'); % Adjust threshold as needed
% disp(['peaks in spectrum ' num2str(periods(peak_locs)) ' sec'])
% %plot it
% figure(3)
% subplot(2,1,1)
% plot(out.t,out.p_c)
% xlabel('Time (s)', 'FontSize', 16)
% ylabel('p_c', 'FontSize', 16)
% subplot(2,1,2)
% findpeaks(spectrum,Fv,'Annotate','extents')
% xlabel('freq (Hz)', 'FontSize', 16)
% ylabel('amplitude', 'FontSize', 16)
% %semilogx(Fv, abs(spectrum),Fv, detrend(abs(spectrum)))
% set(gca,'Xscale','log')
% hold on

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
% disp(['Organ Pipe open-open 1st mode based on mean c is ' num2str(OrganPipe_OO) ' sec'])



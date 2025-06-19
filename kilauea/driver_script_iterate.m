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

% lengths = [700,750,800,850,900,950,1000];
% radi = [5, 10, 15]; 
% n_tots = [0.001, 0.002, 0.003, 0.004, 0.005];
% H2Ofracs = [0.2, 0.4, 0.6, 0.8, 1.0];
% 
lengths = [400, 500, 600, 700];
radi = [10]; 
n_tots = [0.002];
H2Ofracs = [0.4];



loop_params = struct();

for L=1:length(lengths)

    loop_params.L = lengths(L);
    
    for NT=1:length(n_tots)
        
        loop_params.n_tot = n_tots(NT);

        for H=1:length(H2Ofracs)
            
            loop_params.H2Ofrac = H2Ofracs(H);

            for r=1:length(radi)
                loop_params.R = radi(r);
            
                %call script to especify conduit parameters, based on BGstate
                [Mc] = setparameters(loop_params);
                %build the model 
                Model = conduit_internal_g(Mc);
                
                %% time domain simulation.
                
                skip = 1; %only save output every "skip" steps to save memory
                
                T = 1000; %total time in sec
                
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
                
                fun = @(u,t) A*u + Model.Fp(:,1)*Model.M.G(t);
                tic
                
                % Storage arrays
                %ICs
                out.p(:,1) = zeros((Mc.nz+1),1);
                out.t(1) = 0; 
                out.z = Mc.z;
                out.dt = dt;
                out.skip = skip;
                
                out.M = Mc;
                
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
                
                % plotsolutionfields(out)
                
                
                % Make T and Q estimations based on model output
                [T_condres, Q_condres, T_ac, Q_ac] = post_process(out.t,out.p_c,out.dt);
                
                % Identify peaks in spectrum
                %[peaks, peak_locs]=findpeaks(spectrum,'MinPeakProminence',1.5);
                %[peaks, peak_locs] = findpeaks(spectrum2, 'MinPeakHeight', 5); % Adjust threshold as needed
                % [peaks, peak_locs] = findpeaks(spectrum, 'SortStr','descend'); % Adjust threshold as needed
                % disp(['peaks in spectrum ' num2str(periods(peak_locs)) ' sec'])
                
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

                % Extract required variables
                T_est_condres = T_condres;
                T_est_ac = T_ac;
                Q_est_condres = Q_condres;
                Q_est_ac = Q_ac;
                T_RO_model = CRout.T;
                Q_RO_model = CRout.Q;
                R = out.M.R;
                L = out.M.L;
                T_C = out.M.T;
                delta_rho = out.M.rho(1) - out.M.rho(end);
                rho_bar = out.M.rhobar;
                rho_mean = mean(out.M.rho);
                c_mean = mean(out.M.c);
                mu_mean = mean(out.M.mu);
                ngas_mean = mean(out.M.n_gas);
                n_tot = out.M.n_tot;
                H2Ofrac = out.M.H2Ofrac;
                T_ac_analytic = OrganPipe_OO;

                % Define CSV filename
                csv_filename = 'model_results_NEW.csv';
                
                % Check if the file exists; if not, write headers
                if ~isfile(csv_filename)
                    fid = fopen(csv_filename, 'w'); % 'w' mode to create the file and write headers
                    fprintf(fid, 'H2Ofrac,length,radius,n_tot,T_est_condres,T_est_ac,Q_est_condres,Q_est_ac,T_RO_model,Q_RO_model,R,L,T_C,delta_rho,rho_bar,rho_mean,c_mean,mu_mean,ngas_mean,T_ac_analytic\n');
                    fclose(fid);
                end
                
                % Open file in append mode
                fid = fopen(csv_filename, 'a');
                fprintf(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
                    H2Ofrac, loop_params.L, loop_params.R, loop_params.n_tot, ...
                    T_est_condres, T_est_ac, Q_est_condres, Q_est_ac, ...
                    T_RO_model, Q_RO_model, R, L, T_C, delta_rho, ...
                    rho_bar, rho_mean, c_mean, mu_mean, ngas_mean, T_ac_analytic);
                fclose(fid);
                
                % Print a confirmation message
                fprintf('Run finished: H2Ofrac=%.2f, Length=%d, Radius=%d, n_tot=%.3f\n', ...
                    H2Ofrac, loop_params.L, loop_params.R, loop_params.n_tot);

            end
        end
    end
end
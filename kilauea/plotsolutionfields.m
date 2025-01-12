function plotsolutionfields(out)

% plot
grey = [.95 .95 .95];
red = "#BF4539";
orange = "#D9863D";
green = "#BCBF65";
blue = "#51A6A6";
pink = "#BF7582";

figure(1);
subplot(4,1,1);
plot(out.t,out.h,'color', orange, 'linewidth', 1.5);xlabel('time');ylabel('lava lake height (m)')
set(gca,'color',grey)
grid
subplot(4,1,2);
plot(out.t,out.p(round(end/2),:), 'color', blue, 'linewidth', 1.5);xlabel('time');ylabel('conduit pressure (Pa)')
%hold on
%plot(out.t,6*out.p(2,:), 'b', 'linewidth', 1.5);xlabel('time');%ylabel(['column pressure at z=' num2str(round(out.z(round(end/2)))) ' m'])
set(gca,'color',grey)
legend('halfway, L/2')%,'conduit bottom x 6')
grid
if strcmp(out.M.BCtype,'quasistatic')
subplot(4,1,3)
plot(out.t,out.p_c, 'color', pink, 'linewidth', 1.5);xlabel('time');ylabel('P_c (Pa)')
set(gca,'color',grey)
grid
hold on
end

subplot(4,1,4)
hold off
semilogx(out.periods*out.skip, out.spectrum, 'color', blue, 'linewidth', 1)
% if strcmp(out.M.BCtype,'quasistatic')
hold on
semilogx(out.periods*out.skip, out.spectrum2, 'color', red, 'linewidth', 1)
% end
set(gca,'color',grey)
xlim([.5 100])
grid on
xlabel('Period (s)')
ylabel('amplitude spectrum')
legend('Chamber pressure','Conduit halfway down')

% figure(2)
% subplot(1,3,1)
% plot(out.M.R,out.z)
% xlabel('conduit radius (m)')
% ylabel('depth (m)')
% subplot(1,3,2)
% plot(out.M.rho,out.z)
% xlabel('density (kg/m3)')
% title('Background state')
% subplot(1,3,3)
% plot(out.M.c,out.z)
% xlabel('sound speed (m/s)')

figure(3)
pcolor(out.t,out.z,out.p);shading flat
xlabel('time (s)'); ylabel('depth')
title('pressure (Pa)')
cmap
caxis([-1e4,1e4]); colorbar

%some other plots from previous scrips


%xlabel('T (s)')
%ylabel('Amplitude of P_c')
%set(gcf,'color','w')
% hold on;
% for i = 1:length(peak_locs)
%     text(periods(peak_locs(i)), spectrum(peak_locs(i)), sprintf('%.2f seconds', periods(peak_locs(i))), ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     plot(periods(peak_locs(i)), spectrum(peak_locs(i)), 'ro', 'MarkerFaceColor', 'r'); % Mark peak with red circle
% end
% hold off;

% % scalogram
% 
% fb = cwtfilterbank('SignalLength',L,'SamplingPeriod',seconds(1/Fs), 'PeriodLimits', [seconds(1),seconds(50)]);
% scaleSpectrum(fb,out.p_c)
% 


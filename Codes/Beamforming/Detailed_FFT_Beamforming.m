%% =====================================================================
%  FFT-DOMAIN (TRUE-TIME-DELAY) BEAMFORMING  --  step-by-step version
%  ---------------------------------------------------------------------
%  Pipeline:
%    1. Take a wideband signal s(t) and compute its COMPLETE FFT S(f)
%       (magnitude spectrum + phase spectrum).
%    2. N-element array -> which directions can we steer to?
%         * ANY direction (-90..90 deg): the delays can take any value.
%         * N of them are INDEPENDENT (non-overlapping, DFT) beams.
%    3. Direction -> delay:  EVERY ANTENNA gets its own delay
%           tau_n(theta) = n*d*sin(theta)/c ,  n = 0..N-1
%       -> for the N independent beams: N directions x N antennas table.
%    4. Simulate the array:  X_n(f) = S(f) * exp(-j*2*pi*f*tau_n)
%    5. Beamform in the FFT domain (delay = linear phase):
%           W_n(f) = (1/N) * exp(+j*2*pi*f*tau_n(theta_s))
%           Y(f)   = sum_n W_n(f) * X_n(f)   (no inverse FFT is used)
%    6. Change the steering delay, see where output power peaks -> DoA.
%    7. Frequency response of the beamformer, true-delay vs phase-shift.
%    8. Repeat for random target angles.
%  No Phased Array Toolbox needed. Requires MATLAB R2018b+ (xline).
%
%  FORMULAS USED (d = lambda/2, theta measured from broadside)
%    delay at antenna n (ANY theta)  tau_n     = n*d*sin(theta)/c
%    phase step per antenna ........ psi       = 2*pi*f*d*sin(theta)/c  (= pi*sin(theta) at f0)
%    N independent (DFT) beams ..... psi_k     = 2*pi*k/N
%    their directions .............. sin(th_k) = 2k/N ,  k = -N/2 ... N/2-1
%    delay STEP of beam k .......... dtau_k    = d*sin(th_k)/c = k/(N*f0)
%                                    (difference between NEIGHBOURING antennas)
%    delay of antenna n, beam k .... tau_n,k   = n*dtau_k     (8 x 8 table, STEP 3)
%    array input ................... X_n(f)    = S(f)*exp(-j*2*pi*f*tau_n)
%    beamformer .................... Y(f)      = sum_n (1/N)*exp(+j*2*pi*f*tau_steer,n) * X_n(f)
%
%  NOTE on "delay" vs "advance": the weight exp(+j2 pi f tau) is a time
%  ADVANCE. In hardware each antenna instead gets the delay (T - tau_n),
%  with T a common delay chosen so all delays are >= 0. A common delay on
%  every antenna does not change the beam direction, so both are the same
%  beamformer (the table in STEP 3 shows the non-negative version).
%
%  NOTE: k = -4 is theta = -90 deg = END-FIRE (along the array), not usable
%  in practice -> 8 independent beams, 7 usable. Directions between the
%  beams are allowed, but those beams overlap the main 8.
%  =====================================================================
clear; clc; close all;
rng('shuffle');
fprintf('=========== FFT BEAMFORMING: STEP BY STEP ===========\n\n');

%% ---------------------------------------------------------------------
%% STEP 0 : PARAMETERS
%% ---------------------------------------------------------------------
c      = 299792458;          % speed of light (m/s)
f0     = 1.5e9;              % centre frequency (Hz), band = 1-2 GHz
lambda = c/f0;               % wavelength (m)
N      = 8;                  % number of array elements
d      = lambda/2;           % element spacing (m)
pos    = (0:N-1).' * d;      % element positions along the array (N x 1)

fs     = 8e9;                % sampling frequency (Hz)
Ns     = 2048;               % number of samples (= FFT length)
t      = (0:Ns-1).' / fs;    % time vector (s)
df     = fs/Ns;              % FFT bin spacing (Hz)
f      = (-Ns/2:Ns/2-1).' * df;   % two-sided frequency axis (Hz), FFT-shifted
band   = abs(f) >= 1e9 & abs(f) <= 2e9;   % signal band (both sides), used for power

trueAngle  = 25;             % actual direction of the target (deg)
steerAngle = 25;             % beam steering direction for the demo (deg)
steerOff   = 33;             % a deliberately mis-steered beam (deg)
noiseStd   = 0;              % set e.g. 0.05 to add sensor noise

fprintf('f0 = %.2f GHz, lambda = %.2f cm, N = %d, d = lambda/2 = %.2f cm\n', ...
    f0/1e9, lambda*100, N, d*100);
fprintf('fs = %.1f GHz, Ns = %d samples, bin spacing = %.3f MHz\n', ...
    fs/1e9, Ns, df/1e6);
fprintf('Note: at 2 GHz, d = %.3f lambda -> grating lobes if steering beyond about %.0f deg.\n\n', ...
    d*2e9/c, asind(c/(2e9*d) - 1));

%% ---------------------------------------------------------------------
%% STEP 1 : TRANSMIT SIGNAL  s(t)  AND ITS COMPLETE FFT  S(f)
%% ---------------------------------------------------------------------
t0    = t(Ns/2+1);                    % pulse centre (sample Ns/2+1)
sigma = 0.8e-9;                       % envelope width (s)
s     = exp(-0.5*((t-t0)/sigma).^2) .* cos(2*pi*f0*(t-t0));   % Ns x 1, real

% ifftshift moves the pulse centre to sample 1 (t = 0), so there is no
% e^{-j2 pi f t0} term and the phase spectrum is flat. (ifftshift is only a
% circular re-ordering of the samples, NOT an inverse FFT.)
S    = fftshift(fft(ifftshift(s)));
magS = abs(S);
S_ph = complex(real(S), imag(S).*(abs(imag(S)) > 1e-9*max(magS)));  % drop round-off
phS  = angle(S_ph);
phS(magS < 1e-3*max(magS)) = NaN;     % phase is meaningless where |S|~0

[~, k0] = min(abs(f - f0));           % FFT bin closest to f0
fprintf('STEP 1: FFT of s(t).  At f = %.4f GHz (bin %d):\n', f(k0)/1e9, k0);
fprintf('        S(f0) = %.4f %+.4fj   |S| = %.4f   phase = %.2f deg\n\n', ...
    real(S(k0)), imag(S(k0)), abs(S(k0)), rad2deg(angle(S(k0))));

[pk, ipk] = max(magS .* (f > 0));     % spectral peak (positive side)

figure('Name','Step 1: signal and its full FFT','Position',[30 50 1300 700]);
subplot(2,3,1);
plot((t-t0)*1e9, s, 'b'); grid on;
xlim([-8 8]); xlabel('time (ns)'); ylabel('s(t)');
title('(a) Transmit pulse s(t), zoomed around pulse');
subplot(2,3,2);
plot(f/1e9, magS, 'b', 'LineWidth', 1.2); hold on; grid on;
plot(f(ipk)/1e9, pk, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(-f0/1e9,'k:'); xline(f0/1e9,'k:');
xlim([-3 3]); xlabel('frequency (GHz)'); ylabel('|S(f)|');
title('(b) Whole spectrum |S(f)| (two-sided)');
legend('|S(f)|', sprintf('peak %.3f GHz', f(ipk)/1e9), 'Location','north');
subplot(2,3,3);
zz = abs(f - f0) <= 40e6;
stem(f(zz)/1e6, magS(zz), 'b', 'filled', 'MarkerSize', 4); hold on; grid on;
plot(f(ipk)/1e6, pk, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('frequency (MHz)'); ylabel('|S(f)|');
title(sprintf('(c) Zoom: FFT bins around f_0 (spacing %.2f MHz)', df/1e6));
subplot(2,3,[4 6]);
plot(f/1e9, rad2deg(phS), 'r.', 'MarkerSize', 5); grid on;
xlim([-3 3]); ylim([-200 200]); yticks(-180:90:180);
xlabel('frequency (GHz)'); ylabel('\angle S(f) (deg)');
title('(d) Phase spectrum \angle S(f)  (masked where |S| is negligible; pulse centred at t = 0 -> flat at 0^\circ)');

%% ---------------------------------------------------------------------
%% STEP 2 : WHICH DIRECTIONS CAN AN N-ELEMENT ARRAY STEER TO?
%% ---------------------------------------------------------------------
% ANY direction: the delay step d*sin(theta)/c can take any value.
% But only N beams are INDEPENDENT (orthogonal, non-overlapping):
%       psi = pi*sin(theta),  psi_k = 2*pi*k/N  ->  sin(theta_k) = 2k/N
% Each of these has a null in the directions of the other N-1 beams.
kk         = (-N/2:N/2-1).';
sinTheta_k = 2*kk/N;
theta_k    = asind(sinTheta_k);
dtau_k     = d*sinTheta_k/c;          % delay STEP between neighbouring antennas

fprintf('STEP 2: the beam can point in ANY direction from -90 to +90 deg.\n');
fprintf('        N = %d antennas give N = %d INDEPENDENT (DFT) beams:\n', N, N);
fprintf('   k   sin(theta_k)   theta_k(deg)   phase step psi_k(deg)   delay STEP between neighbours (ps)\n');
for i = 1:N
    fprintf('  %2d     %+6.3f        %+7.2f            %+7.1f                   %+7.2f\n', ...
        kk(i), sinTheta_k(i), theta_k(i), rad2deg(2*pi*kk(i)/N), dtau_k(i)*1e12);
end
fprintf('   The delay STEP is not the delay of one antenna: antenna n gets n x step (see STEP 3).\n');
fprintf('   Beams are uniform in sin(theta), not theta -> wider near +/-90 deg.\n');
fprintf('   Directions in between (e.g. %d deg) are allowed; they overlap the main beams.\n', trueAngle);
fprintf('   k = -4 (-90 deg) is END-FIRE (along the array): 8 independent beams, 7 usable.\n\n');

figure('Name','Q1: steering directions','Position',[60 60 1100 480]);
subplot(1,2,1); hold on;
plot([-0.9 0.9], [-0.08 -0.08], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
plot(linspace(-0.9,0.9,N), -0.08*ones(1,N), 'ks', 'MarkerFaceColor','k', 'MarkerSize', 5);
text(0, -0.2, 'array (8 elements)', 'HorizontalAlignment','center');
for i = 1:N
    if abs(theta_k(i)) > 89.9
        clr = [0.55 0.55 0.55];  lab = sprintf('k=%d\n%.1f^\\circ\n(end-fire)', kk(i), theta_k(i));
    else
        clr = [0 0.3 0.9];       lab = sprintf('k=%d\n%.1f^\\circ', kk(i), theta_k(i));
    end
    plot([0 sind(theta_k(i))], [0 cosd(theta_k(i))], '-', 'Color', clr, 'LineWidth', 2.2);
    text(1.18*sind(theta_k(i)), 1.18*cosd(theta_k(i)), lab, ...
        'HorizontalAlignment','center', 'FontSize', 9, 'Color', clr);
end
plot([0 sind(trueAngle)], [0 cosd(trueAngle)], '--', 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
text(1.05*sind(trueAngle)+0.25, 1.05*cosd(trueAngle), sprintf('%d^\\circ (in between)', trueAngle), ...
    'Color', [0.8 0.2 0.2], 'FontSize', 9);
axis equal; axis([-1.6 1.6 -0.3 1.4]); axis off;
title({'8 independent beams (blue = 7 usable, grey = end-fire)', ...
       'red dashed = any other direction is also possible'});
subplot(1,2,2);
thC = -90:0.5:90;
plot(thC, d*sind(thC)/c*1e12, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); hold on; grid on;
plot(theta_k, dtau_k*1e12, 'o', 'Color', [0 0.3 0.9], 'MarkerFaceColor', [0 0.3 0.9]);
plot(trueAngle, d*sind(trueAngle)/c*1e12, 's', 'Color', [0.8 0.2 0.2], 'MarkerFaceColor', [0.8 0.2 0.2]);
xlabel('direction \theta (deg)'); ylabel('delay step between neighbours (ps)');
title('Delay step = d sin\theta / c  (can take ANY value)');
legend('any direction','8 independent beams', sprintf('example %d^\\circ', trueAngle), 'Location','northwest');

%% ---------------------------------------------------------------------
%% STEP 3 : DIRECTION -> DELAY OF EVERY ANTENNA  tau_n = n*d*sin(theta)/c
%% ---------------------------------------------------------------------
tau_true  = pos*sind(trueAngle)/c;        % N x 1  (true arrival delays)
tau_steer = pos*sind(steerAngle)/c;       % N x 1  (steering delays)

fprintf('STEP 3: EACH ANTENNA has its own delay.  tau_n = n*d*sin(theta)/c\n');
fprintf('   Example, theta = %d deg (not one of the 8 beams):\n', trueAngle);
fprintf('   n    pos (lambda)   tau_true (ps)   tau_steer (ps)\n');
for n = 1:N
    fprintf('  %2d      %5.2f         %8.2f        %8.2f\n', ...
        n-1, pos(n)/lambda, tau_true(n)*1e12, tau_steer(n)*1e12);
end
fprintf('\n');

figure('Name','Q2: delay of each antenna','Position',[60 60 700 450]);
stem(0:N-1, tau_true*1e12, 'filled'); grid on;
xlabel('antenna number n'); ylabel('delay \tau_n (ps)');
title(sprintf('Delay of each antenna, target at %d^\\circ  (\\tau_n = n d sin\\theta / c)', trueAngle));

% ---- 8 directions x 8 antennas table ----
% Row i = direction theta_k(i), column n = antenna n-1:   tau = n * dtau_k
tauTable = dtau_k .* (0:N-1);                          % N x N (s), may be negative
tauPhys  = tauTable - min(tauTable, [], 2);            % add a common delay per row so all >= 0
                                                       % (same shift on all antennas -> same beam)
fprintf('Delay of EVERY antenna for each of the %d independent directions (ps):\n', N);
fprintf(' %3s   %7s   |', 'k', 'theta');
for n = 0:N-1, fprintf(' %7s', sprintf('ant%d', n)); end
fprintf('\n');
for i = 1:N
    fprintf(' %3d   %+7.1f   |', kk(i), theta_k(i));
    fprintf(' %7.1f', tauPhys(i,:)*1e12);
    fprintf('\n');
end
tauEx = (0:N-1)*d*sind(trueAngle)/c;
fprintf(' %3s   %+7.1f   |', '-', trueAngle); fprintf(' %7.1f', tauEx*1e12);
fprintf('   <- any other direction: one more row of 8\n\n');

figure('Name','Delay of each antenna, for each direction','Position',[80 80 850 500]);
cl = lines(N);  cl(N,:) = [0.1 0.1 0.1];    % lines(8) repeats colour 1 as colour 8 -> replace it
hold on;
for i = 1:N
    lw = 1.4;  ls = '-o';
    if abs(theta_k(i)) > 89.9, cl(i,:) = [0.55 0.55 0.55]; ls = '--o'; end
    if kk(i) == 0, lw = 2.8; end
    plot(0:N-1, tauPhys(i,:)*1e12, ls, 'Color', cl(i,:), 'LineWidth', lw, 'MarkerSize', 4);
end
plot(0:N-1, tauEx*1e12, ':s', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.6, 'MarkerSize', 4);
grid on; xlabel('antenna number n'); ylabel('delay given to the antenna (ps)');
title({'Each antenna gets its own delay', ...
       'x: antenna number     y: delay (ps)     one line = one direction'});
legend([arrayfun(@(a) sprintf('%.1f^\\circ', a), theta_k, 'UniformOutput', false); ...
        {sprintf('%d^\\circ (example)', trueAngle)}], 'Location','northeastoutside');

%% ---------------------------------------------------------------------
%% STEP 4 : SIMULATE THE ARRAY  X_n(f) = S(f) * exp(-j 2 pi f tau_n)
%% ---------------------------------------------------------------------
X = receiveArray(S, f, pos, c, trueAngle, noiseStd, Ns);   % Ns x N

relPhaseMeas = angle(X(k0,:) ./ X(k0,1)).';
relPhaseTheo = angle(exp(-1j*2*pi*f(k0)*(tau_true - tau_true(1))));
fprintf('STEP 4: received element spectra X_n(f) = S(f) e^{-j2 pi f tau_n}\n');
fprintf('   Phase of X_n(f0) relative to element 0  (= -2 pi f0 tau_n, wrapped):\n');
fprintf('   n    measured (deg)   theory (deg)\n');
for n = 1:N
    fprintf('  %2d      %8.2f        %8.2f\n', n-1, rad2deg(relPhaseMeas(n)), rad2deg(relPhaseTheo(n)));
end
fprintf('\n');

bi = f >= 1e9 & f <= 2e9;
figure('Name','Step 4: array received signals','Position',[80 60 1000 750]);
subplot(2,2,1);
stem(0:N-1, tau_true*1e12, 'filled'); grid on;
xlabel('element n'); ylabel('\tau_n (ps)');
title(sprintf('True arrival delays, \\theta = %d^\\circ', trueAngle));
subplot(2,2,2);
stem(0:N-1, rad2deg(relPhaseMeas), 'filled'); hold on;
plot(0:N-1, rad2deg(relPhaseTheo), 'ro');  grid on;
xlabel('element n'); ylabel('phase (deg)');
title('Phase of X_n(f_0) rel. to element 0: measured (stem) vs theory (o)');
subplot(2,2,3);
plot(f(bi)/1e9, unwrap(angle(X(bi,:) ./ S(bi))), 'LineWidth', 1.1); grid on;
xlabel('frequency (GHz)'); ylabel('phase of X_n/S (rad)');
title('Phase of each element vs f: slope = -2\pi\tau_n');
legend(arrayfun(@(n) sprintf('n=%d',n), 0:N-1, 'UniformOutput', false), 'Location','southwest','NumColumns',2);
subplot(2,2,4);
plot(f/1e9, abs(X)); grid on; xlim([-3 3]);
xlabel('frequency (GHz)'); ylabel('|X_n(f)|');
title('Element magnitude spectra (delay only changes phase)');

%% ---------------------------------------------------------------------
%% STEP 5 : BEAMFORMING  W_n(f) = (1/N) e^{+j 2 pi f tau_steer,n}
%%                       Y(f)   = sum_n W_n(f) X_n(f)
%% ---------------------------------------------------------------------
W = exp(+1j*2*pi*f*tau_steer.') / N;         % Ns x N steering weights (advance, see header)
Y = sum(W .* X, 2);                          % Ns x 1  beamformer output spectrum

fprintf('STEP 5: steering weights at f0 for theta_s = %d deg\n', steerAngle);
fprintf('   w_n = (1/N) exp(+j 2 pi f0 tau_steer,n)\n');
fprintf('   n    |w_n|    phase(w_n) deg    phase(w_n X_n) rel. elem 0 (deg)\n');
prod0 = W(k0,:) .* X(k0,:);
for n = 1:N
    fprintf('  %2d   %.4f     %9.2f            %9.2f\n', n-1, abs(W(k0,n)), ...
        rad2deg(angle(W(k0,n))), rad2deg(angle(prod0(n)/prod0(1))));
end
fprintf('   -> after weighting every element has the SAME phase (last column ~ 0),\n');
fprintf('      so the N terms add coherently:  |Y(f0)| = %.4f   vs   |S(f0)| = %.4f\n\n', ...
    abs(Y(k0)), abs(S(k0)));

figure('Name','Step 5: steering weight vectors','Position',[110 70 950 700]);
subplot(2,2,1);
stem(0:N-1, abs(W(k0,:)), 'filled'); grid on; ylim([0 0.2]);
xlabel('element n'); ylabel('|w_n|'); title('Weight magnitude at f_0 (= 1/N)');
subplot(2,2,2);
stem(0:N-1, rad2deg(angle(W(k0,:))), 'filled'); grid on;
xlabel('element n'); ylabel('phase (deg)');
title(sprintf('Weight phase at f_0 (wrapped), \\theta_s = %d^\\circ', steerAngle));
subplot(2,2,[3 4]);
imagesc(f(bi)/1e9, 0:N-1, rad2deg(angle(W(bi,:))).'); axis xy; colormap(gca,'jet');
cb = colorbar; ylabel(cb,'phase (deg)');
xlabel('frequency (GHz)'); ylabel('element n');
title('Weight phase of every element at every frequency (phase ramp grows with f)');

W_off = exp(+1j*2*pi*f*(pos*sind(steerOff)/c).') / N;
Y_off = sum(W_off .* X, 2);
W_0   = ones(Ns,N)/N;                        % no delays at all (steer 0 deg)
Y_0   = sum(W_0 .* X, 2);

fprintf('Peak |spectrum|: single element = %.3f | steered %d deg (aligned) = %.3f | steered %d deg = %.3f | no delay = %.3f\n\n', ...
    max(abs(X(:,1))), steerAngle, max(abs(Y)), steerOff, max(abs(Y_off)), max(abs(Y_0)));

figure('Name','Step 5: beamformer output','Position',[140 80 1000 750]);
subplot(3,1,1);
plot(f/1e9, abs(S), 'k', 'LineWidth', 1.5); hold on;
plot(f/1e9, abs(X(:,1)), 'g--');
plot(f/1e9, abs(Y), 'r', 'LineWidth', 1.1); grid on; xlim([-3 3]);
legend('|S(f)| original','|X_0(f)| one element','|Y(f)| beamformer','Location','northeast');
xlabel('frequency (GHz)'); ylabel('magnitude'); title('Output spectrum magnitude');
subplot(3,1,2);
mS = abs(S) > 1e-3*max(abs(S));  mY = abs(Y) > 1e-3*max(abs(Y));
plot(f(mS)/1e9, rad2deg(angle(S(mS))), 'k.', 'MarkerSize', 4); hold on;
plot(f(mY)/1e9, rad2deg(angle(Y(mY))), 'r.', 'MarkerSize', 4); grid on; xlim([-3 3]);
legend('\angle S(f)','\angle Y(f)'); xlabel('frequency (GHz)'); ylabel('phase (deg)');
title('Output phase spectrum: identical to S(f) when steering = target direction');
subplot(3,1,3);
plot(f/1e9, abs(S), 'k', 'LineWidth', 1.5); hold on;
plot(f/1e9, abs(Y), 'r', 'LineWidth', 1.1);
plot(f/1e9, abs(Y_off), 'm:', 'LineWidth', 1.3);
plot(f/1e9, abs(Y_0), 'c-.');
grid on; xlim([-3 3]); xlabel('frequency (GHz)'); ylabel('|Y(f)|');
legend('|S(f)| original','steered at target', ...
    sprintf('steered %d^\\circ (miss)', steerOff), 'no delay', 'Location','northeast');
title('Output spectrum |Y(f)| for different steering delays');

%% ---------------------------------------------------------------------
%% STEP 6 : CHANGE THE DELAY, SEE WHICH DIRECTION THE BEAM POINTS AT
%% ---------------------------------------------------------------------
scanStep  = 0.01;                                   % scan step (deg): make smaller for more directions
thetaScan = -90:scanStep:90;                         % ANY direction
fprintf('Scan step = %g deg -> %d steering directions\n', scanStep, numel(thetaScan));
P_scan    = scanPower(X, f, pos, c, thetaScan, band);
[~, im]   = max(P_scan);
estAngle  = thetaScan(im);

P_disc  = scanPower(X, f, pos, c, theta_k.', band);  % only the 8 independent beams
[~, ik] = max(P_disc);

fprintf('STEP 6: sweeping the steering delay over %d directions\n', numel(thetaScan));
fprintf('   Any delay (fine sweep): power peaks at %.2f deg  (true = %d deg)\n', estAngle, trueAngle);
fprintf('   Only the %d fixed beams: best beam k = %d -> %.1f deg\n', N, kk(ik), theta_k(ik));
fprintf('   (fixed-beam grid limits RESOLUTION, not where the array can point)\n\n');

figure('Name','Q3: output vs steering delay','Position',[60 60 1100 480]);
subplot(1,2,1);
plot(thetaScan, P_scan/max(P_scan), 'b', 'LineWidth', 1.5); hold on; grid on;
plot(estAngle, 1, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(trueAngle, 'k--', 'true direction');
xlabel('steering angle (deg)'); ylabel('output power (normalised)');
title(sprintf('Any delay (fine sweep): peak at %.2f^\\circ', estAngle));
subplot(1,2,2);
bar(kk, P_disc/max(P_disc), 'FaceColor', [0.6 0.75 0.95]); hold on;
bar(kk(ik), P_disc(ik)/max(P_disc), 'FaceColor', [0.85 0.2 0.2]);
grid on; xticks(kk); xticklabels(arrayfun(@(a) sprintf('%.1f',a), theta_k, 'UniformOutput', false));
xlabel('beam direction (deg)'); ylabel('output power (normalised)');
title(sprintf('Only the 8 fixed beams: nearest is %.1f^\\circ', theta_k(ik)));

%% ---------------------------------------------------------------------
%% STEP 7 : FREQUENCY RESPONSE OF THE BEAMFORMER
%% ---------------------------------------------------------------------
%   True-time-delay (FFT) beamformer:
%       H(f,theta) = (1/N) sum_n exp( j 2 pi f n d (sin(theta_s)-sin(theta))/c )
%   Single phase-shift beamformer (weights fixed at f0 only):
%       H(f,theta) = (1/N) sum_n exp( j 2 pi n d (f0 sin(theta_s) - f sin(theta))/c )
thetaG = -90:0.25:90;
fPlot  = [1.0 1.5 2.0]*1e9;
Hpat_ttd = zeros(numel(fPlot), numel(thetaG));
Hpat_ps  = zeros(numel(fPlot), numel(thetaG));
for i = 1:numel(fPlot)
    Hpat_ttd(i,:) = abs(mean(exp(1j*2*pi*fPlot(i)/c*(pos*(sind(steerAngle)-sind(thetaG)))), 1));
    Hpat_ps(i,:)  = abs(mean(exp(1j*2*pi/c*(pos*(f0*sind(steerAngle)-fPlot(i)*sind(thetaG)))), 1));
end

fb   = linspace(0.5e9, 2.5e9, 801).';
st_t = sind(trueAngle);
H_ttd_on = mean(exp(1j*2*pi*fb/c*(pos.'*(sind(steerAngle)-st_t))), 2);
H_ps_on  = mean(exp(1j*2*pi/c*((f0*sind(steerAngle)-fb*st_t)*pos.')), 2);
H_ttd_of = mean(exp(1j*2*pi*fb/c*(pos.'*(sind(steerOff)-st_t))), 2);

mk    = f > 0.5e9 & f < 2.5e9 & abs(S) > 0.02*max(abs(S));
Hmeas = Y_off ./ S;

fprintf('STEP 7: frequency response computed (see figure).\n');
fprintf('   True-delay beamformer on target: |H(f)| = 1 at all f (no distortion).\n');
fprintf('   Single phase-shift beamformer:   |H(f)| = %.3f at 1.0 GHz, %.3f at 2.0 GHz\n\n', ...
    abs(interp1(fb, H_ps_on, 1.0e9)), abs(interp1(fb, H_ps_on, 2.0e9)));

figure('Name','Step 7: frequency response','Position',[200 60 1100 800]);
subplot(2,2,1);
plot(thetaG, 20*log10(Hpat_ttd'+eps), 'LineWidth', 1.2); grid on; ylim([-40 2]);
xline(steerAngle,'k--'); xlabel('\theta (deg)'); ylabel('|H| (dB)');
title('True-delay (FFT) beamformer: pattern at 3 frequencies');
legend('1.0 GHz','1.5 GHz','2.0 GHz','Location','south');
subplot(2,2,2);
plot(thetaG, 20*log10(Hpat_ps'+eps), 'LineWidth', 1.2); grid on; ylim([-40 2]);
xline(steerAngle,'k--'); xlabel('\theta (deg)'); ylabel('|H| (dB)');
title('Single phase-shift (f_0 only): beam squints with frequency');
legend('1.0 GHz','1.5 GHz','2.0 GHz','Location','south');
subplot(2,2,3);
plot(fb/1e9, abs(H_ttd_on), 'b', 'LineWidth', 1.5); hold on;
plot(fb/1e9, abs(H_ps_on),  'r', 'LineWidth', 1.5);
plot(fb/1e9, abs(H_ttd_of), 'g', 'LineWidth', 1.5);
plot(f(mk)/1e9, abs(Hmeas(mk)), 'ko', 'MarkerSize', 3);
plot(f/1e9, abs(S)/max(abs(S)), 'Color',[0.6 0.6 0.6]);
grid on; xlim([0.5 2.5]); ylim([0 1.1]);
xlabel('frequency (GHz)'); ylabel('|H(f)|');
title(sprintf('Response for source at \\theta = %d^\\circ', trueAngle));
legend('true-delay, steered on target','phase-shift, steered on target', ...
    sprintf('true-delay, steered %d^\\circ',steerOff),'measured Y/S (simulation)','|S(f)| (scaled)', ...
    'Location','southwest');
subplot(2,2,4);
fg = linspace(0.5e9, 2.5e9, 200);  tg = -90:1:90;
Hmap = zeros(numel(fg), numel(tg));
for i = 1:numel(fg)
    Hmap(i,:) = abs(mean(exp(1j*2*pi/c*(pos*(f0*sind(steerAngle)-fg(i)*sind(tg)))), 1));
end
imagesc(tg, fg/1e9, 20*log10(Hmap+eps), [-40 0]); axis xy; colorbar;
xlabel('\theta (deg)'); ylabel('frequency (GHz)');
title('Phase-shift beamformer |H(f,\theta)| (dB): ridge bends with f');
hold on; xline(steerAngle,'w--');

%% ---------------------------------------------------------------------
%% STEP 8 : REPEAT FOR SEVERAL RANDOM TARGET ANGLES (moving target)
%% ---------------------------------------------------------------------
% A ULA cannot tell theta from 180-theta, so only -90..90 is scanned.
% Targets are drawn from -80..80 deg: near +/-90 (end-fire) the beam is
% very wide and the estimate is poor.
numTrials  = 8;
trueAngles = randi([-80 80], 1, numTrials);
estAngles  = zeros(1, numTrials);
for tr = 1:numTrials
    Xr = receiveArray(S, f, pos, c, trueAngles(tr), noiseStd, Ns);
    Pr = scanPower(Xr, f, pos, c, thetaScan, band);
    [~, ii] = max(Pr);
    estAngles(tr) = thetaScan(ii);
    fprintf('Trial %d: true = %4d deg, estimated = %6.1f deg, error = %.1f deg\n', ...
        tr, trueAngles(tr), estAngles(tr), abs(trueAngles(tr)-estAngles(tr)));
end
fprintf('\nMean absolute error over %d trials: %.2f deg\n', numTrials, mean(abs(trueAngles-estAngles)));

figure('Name','Step 8: tracking accuracy');
plot(1:numTrials, trueAngles, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(1:numTrials, estAngles,  'rx--', 'LineWidth', 1.5, 'MarkerSize', 10); grid on;
xlabel('trial'); ylabel('angle (deg)'); legend('true','estimated','Location','best');
title('True vs estimated direction (FFT-domain delay-and-sum scan)');

%% =====================  LOCAL FUNCTIONS  ==============================
function X = receiveArray(S, f, pos, c, ang, noiseStd, Ns)
    % Received element spectra for a plane wave from 'ang' degrees
    tau = pos*sind(ang)/c;
    N   = numel(pos);
    X   = S .* exp(-1j*2*pi*f*tau.');
    X   = X + noiseStd*sqrt(Ns/2)*(randn(Ns,N) + 1j*randn(Ns,N));
end

function P = scanPower(X, f, pos, c, grid, band)
    % Output power of the FFT-domain delay-and-sum beamformer for every
    % steering angle in 'grid', summed only over the signal band
    % (out-of-band bins hold only noise).
    N  = numel(pos);
    P  = zeros(size(grid));
    fb = f(band);  Xb = X(band, :);          % keep only the 1-2 GHz bins
    for i = 1:numel(grid)
        tau_s = pos*sind(grid(i))/c;
        Wg    = exp(1j*2*pi*fb*tau_s.') / N;   % only the band bins -> faster
        Yg    = sum(Wg .* Xb, 2);
        P(i)  = sum(abs(Yg).^2);
    end
end
%% STEERING IS NOT LIMITED: many beams on polar plots
%  Runs on its own, or paste at the end of fft_beamforming.m.
if ~exist('N','var'),  N  = 8;          end
if ~exist('c','var'),  c  = 299792458;  end
if ~exist('f0','var'), f0 = 1.5e9;      end
if ~exist('d','var'),  d  = c/f0/2;     end

n    = (0:N-1).';
th   = -90:0.25:90;
beam = @(a) abs(mean(exp(1j*2*pi*f0*n*d*(sind(a)-sind(th))/c), 1));

figure('Name','Steering is not limited','Position',[50 80 1400 520]);

% (1) beam steered every 10 deg from -60 to +60
subplot(1,2,1);
angs = -60:10:60;
cols = turbo(numel(angs));
for i = 1:numel(angs)
    polarplot(deg2rad(th), beam(angs(i)), 'Color', cols(i,:), 'LineWidth', 1.6); hold on;
end
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; thetalim([-90 90]);
legend(arrayfun(@(a) sprintf('%d^\\circ', a), angs, 'UniformOutput', false), ...
    'Location','southoutside','NumColumns',7);
title({'(1) Beam steered every 10^\circ from -60^\circ to +60^\circ', ...
       '13 beams, each points exactly where we set it'});

% (2) zoom: beams IN BETWEEN two of the 8 beams (14.5 and 30 deg)
subplot(1,2,2);
h1 = polarplot(deg2rad(th), beam(asind(0.25)), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 2.2); hold on;
polarplot(deg2rad(th), beam(30), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 2.2);
angs2 = 16:2:28;
cols2 = autumn(numel(angs2)+2);
hh = gobjects(1, numel(angs2));
for i = 1:numel(angs2)
    hh(i) = polarplot(deg2rad(th), beam(angs2(i)), 'Color', cols2(i,:), 'LineWidth', 1.6);
end
ax = gca; ax.ThetaZeroLocation = 'top'; ax.ThetaDir = 'clockwise'; thetalim([-90 90]);
legend([hh h1], [arrayfun(@(a) sprintf('%d^\\circ', a), angs2, 'UniformOutput', false), ...
    {'14.5^\circ & 30^\circ (two of the 8 beams)'}], 'Location','southoutside','NumColumns',4);
title({'(2) Zoom: 7 beams IN BETWEEN 14.5^\circ and 30^\circ', ...
       'steering is not limited to the 8 beams'});
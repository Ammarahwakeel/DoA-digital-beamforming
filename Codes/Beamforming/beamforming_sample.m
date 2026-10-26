%% =====================================================================
%  FFT-BASED DELAY-AND-SUM BEAMFORMING WITH WHOLE-SAMPLE DELAYS
%  ---------------------------------------------------------------------
%  In the digital domain a delay is a shift by a WHOLE number of samples,
%  so every delay must be a multiple of Ts = 1/fs. This limits the beams:
%
%    delay step between neighbouring antennas = m * Ts      (m = integer)
%    beam direction                            sin(theta_m) = m*c*Ts / d
%    |sin(theta)| <= 1   ->   |m| <= floor(d / (c*Ts))
%
%  With fs = 8 GS/s (Ts = 125 ps) and d = 9.99 cm:  |m| <= 2  ->  5 beams
%  (0, +/-22.02, +/-48.59 deg).
%
%  This script
%    1. finds the possible beams and the whole-sample delay of every antenna
%    2. simulates the 8 received time signals for a target at trueAngle
%    3. takes the FFT of every antenna, applies each whole-sample delay as
%       a phase  exp(-j*2*pi*k*a/Ns), adds the antennas (FFT-based
%       delay-and-sum), and checks the result against a plain time shift
%    4. picks the beam with the most output power
%  and writes three files used by the C program (same folder):
%    fftbf_delay_table.csv   delays measured here, in samples
%    fftbf_rx_signals.csv    received time samples, Ns rows x 8 antennas
%    fftbf_matlab_power.csv  output power of every beam (to compare with C)
%  =====================================================================
clear; clc; close all;

%% ---------------- PARAMETERS ------------------------------------------
c      = 299792458;          % speed of light (m/s)
f0     = 1.5e9;              % centre frequency (Hz), band 1-2 GHz
lambda = c/f0;
N      = 8;                  % antennas
d      = lambda/2;           % spacing (m)
fs     = 8e9;                % sampling frequency (Hz)
Ts     = 1/fs;               % sample period (s)
Ns     = 2048;               % samples
t      = (0:Ns-1).'/fs;
trueAngle = 25;              % target direction (deg), unknown to the receiver
n      = 0:N-1;              % antenna index

fprintf('fs = %.1f GS/s -> Ts = %.1f ps, d = %.2f cm, c*Ts = %.2f cm\n', ...
    fs/1e9, Ts*1e12, d*100, c*Ts*100);

%% ---------------- STEP 1: POSSIBLE BEAMS AND THEIR DELAYS -------------
mMax    = floor(d/(c*Ts));             % largest whole-sample step between neighbours
m       = (-mMax:mMax).';              % allowed steps (samples)
thBeam  = asind(m*c*Ts/d);             % beam directions (deg)
nB      = numel(m);
arrival = m * n;                       % nB x N: arrival delay of antenna n (samples)
align   = max(arrival, [], 2) - arrival;   % delay we APPLY to line them up (>= 0)

fprintf('\nSTEP 1: delay step must be a whole number of samples -> |m| <= %d -> %d beams\n', mMax, nB);
fprintf('   m   step (ps)   beam (deg) | delay applied to antenna 0..%d (samples)\n', N-1);
for b = 1:nB
    fprintf('  %+d   %+6.0f     %+7.2f   |', m(b), m(b)*Ts*1e12, thBeam(b));
    fprintf(' %3d', align(b,:));
    fprintf('\n');
end

fid = fopen('fftbf_delay_table.csv', 'w');
fprintf(fid, 'm,theta_deg'); fprintf(fid, ',ant%d', n); fprintf(fid, '\n');
for b = 1:nB
    fprintf(fid, '%d,%.6f', m(b), thBeam(b)); fprintf(fid, ',%d', align(b,:)); fprintf(fid, '\n');
end
fclose(fid);

%% ---------------- STEP 2: RECEIVED SIGNALS -----------------------------
% Transmit pulse (same as the FFT script) and the TRUE arrival delays.
% The true delays are physical, so they are NOT whole samples; they are
% applied exactly (in the frequency domain) only to SIMULATE the antennas.
t0    = t(Ns/2+1);
sigma = 0.8e-9;
s     = exp(-0.5*((t-t0)/sigma).^2) .* cos(2*pi*f0*(t-t0));
tau   = n.' * d * sind(trueAngle) / c;            % N x 1 true delays (s)
fk    = (0:Ns-1).' * fs/Ns;
fk(fk >= fs/2) = fk(fk >= fs/2) - fs;              % signed FFT frequencies
x     = real(ifft(fft(s) .* exp(-1j*2*pi*fk*tau.')));   % Ns x N received signals

fprintf('\nSTEP 2: target at %g deg -> true delays (in samples):', trueAngle);
fprintf(' %.2f', tau/Ts); fprintf('\n');
fprintf('        (not whole numbers: a real target is rarely exactly on a beam)\n');

fid = fopen('fftbf_rx_signals.csv', 'w');
for k = 1:Ns
    fprintf(fid, '%.17g', x(k,1)); fprintf(fid, ',%.17g', x(k,2:end)); fprintf(fid, '\n');
end
fclose(fid);

%% ---------------- STEP 3: FFT-BASED DELAY-AND-SUM ----------------------
% The whole-sample delays are applied in the FFT domain. A delay of a
% whole samples is a phase:   x[t - a]  <->  X(k) * exp(-j*2*pi*k*a/Ns)
% so only the delay a (from the table) and the bin number k are needed.
Xf = fft(x);                                   % Ns x N, FFT of every antenna
kb = (0:Ns-1).';                               % FFT bin index
P  = zeros(nB,1);
Yf = zeros(Ns,nB);
for b = 1:nB
    Wf      = exp(-1j*2*pi*kb*align(b,:)/Ns);  % Ns x N: phase of each antenna's delay
    Yf(:,b) = sum(Wf .* Xf, 2) / N;            % delay-and-sum in the FFT domain
    P(b)    = sum(abs(Yf(:,b)).^2) / Ns;       % output power (Parseval: = time-domain power)
end
[~, ib]  = max(P);
estAngle = thBeam(ib);

% Check: the same delays applied as a time shift give (practically) the same powers
Ptime = zeros(nB,1);
for b = 1:nB, Ptime(b) = sum(delaySum(x, align(b,:)).^2); end
fprintf('\nSTEP 3: FFT-based delay-and-sum (time-shift check: largest difference %.1e of the peak)\n', ...
    max(abs(P - Ptime))/max(P));
fprintf('   beam (deg)    power      normalised\n');
for b = 1:nB
    fprintf('   %+7.2f    %10.4f    %.4f\n', thBeam(b), P(b), P(b)/max(P));
end
fprintf('\nEstimated direction: %.2f deg  (true %g deg, error %.2f deg)\n', ...
    estAngle, trueAngle, abs(estAngle - trueAngle));
fprintf('With whole-sample delays the answer can only be one of the %d beams.\n', nB);

fid = fopen('fftbf_matlab_power.csv', 'w');
fprintf(fid, 'm,theta_deg,power\n');
for b = 1:nB, fprintf(fid, '%d,%.6f,%.17g\n', m(b), thBeam(b), P(b)); end
fclose(fid);
fprintf('Files written: fftbf_delay_table.csv, fftbf_rx_signals.csv, fftbf_matlab_power.csv\n');

%% ---------------- FIGURES ----------------------------------------------
figure('Name','Beam powers','Position',[60 60 800 450]);
bar(thBeam, P/max(P), 0.5); hold on; grid on;
bar(thBeam(ib), 1, 0.5, 'FaceColor', [0.85 0.2 0.2]);
xline(trueAngle, 'k--', sprintf('true %g^\\circ', trueAngle), 'LineWidth', 1.5);
xlim([-90 90]); xlabel('beam direction (deg)'); ylabel('output power (normalised)');
title(sprintf('Whole-sample delays: only %d beams, best = %.2f^\\circ', nB, estAngle));

figure('Name','The possible beams');
th  = -90:0.25:90;
pax = polaraxes; hold(pax, 'on');
for b = 1:nB
    H = abs(mean(exp(1j*2*pi*f0*(n.'*m(b)*Ts - n.'*d*sind(th)/c)), 1));
    polarplot(pax, deg2rad(th), H, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('m = %+d  (%+.2f^\\circ)', m(b), thBeam(b)));
end
pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise'; thetalim(pax, [-90 90]);
legend(pax, 'Location', 'southoutside', 'NumColumns', 3);
title(pax, sprintf('All %d beams possible with whole-sample delays (fs = %g GS/s, 1.5 GHz)', nB, fs/1e9));

figure('Name','Alignment','Position',[80 80 1000 500]);
win = (Ns/2-40:Ns/2+60);
subplot(1,2,1); plot((win-1)/fs*1e9, x(win,:)); grid on;
xlabel('time (ns)'); ylabel('amplitude'); title('Received signals (arrive at different times)');
subplot(1,2,2); hold on; grid on;
for k = 1:N
    a  = align(ib,k);
    xd = [zeros(a,1); x(1:Ns-a,k)];            % antenna k delayed by a samples
    plot((win-1)/fs*1e9, xd(win));
end
xlabel('time (ns)'); ylabel('amplitude');
title(sprintf('After whole-sample delays of the best beam (%.2f^\\circ)', estAngle));

%% ---------------- LOCAL FUNCTION ---------------------------------------
function y = delaySum(x, a)
    % Delay antenna n by a(n) WHOLE samples, add all antennas, divide by N
    [Ns, N] = size(x);
    y = zeros(Ns, 1);
    for k = 1:N
        y(a(k)+1:Ns) = y(a(k)+1:Ns) + x(1:Ns-a(k), k);
    end
    y = y / N;
end
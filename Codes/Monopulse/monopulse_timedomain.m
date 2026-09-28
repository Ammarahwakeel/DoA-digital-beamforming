%% =======================================================================
%  MONOPULSE ON AN 8x1 ULA USING MATLAB BUILT-IN FUNCTIONS
%  Built-ins used:
%     phased.ULA                            -> 8x1 antenna array
%     collectPlaneWave                      -> signal received by each antenna
%     pattern                               -> beam pattern (main / side / back lobes)
%     phased.SumDifferenceMonopulseTracker  -> monopulse angle estimation
%  Pure monopulse only (beams at boresight, no beamforming / steering).
%  =======================================================================
clear; clc; close all;
rng(1);

%% 1. System setup (8x1 array, 1-2 GHz radio band) -----------------------
carrierFreq = 1.5e9;                              % centre of 1-2 GHz band
speed       = physconst('LightSpeed');
wavelength  = speed/carrierFreq;                  % 20 cm
spacing     = (speed/2e9)/2;                      % 7.5 cm (lambda/2 at 2 GHz)
numElements = 8;

ula = phased.ULA('NumElements', numElements, 'ElementSpacing', spacing);

% Built-in monopulse: splits the ULA into two halves, forms Sum and
% Difference channels and converts their ratio into an azimuth angle.
tracker = phased.SumDifferenceMonopulseTracker('SensorArray', ula, ...
    'PropagationSpeed', speed, 'OperatingFrequency', carrierFreq);
steerAng = 0;                                     % boresight -> pure monopulse

nullAngle = asind(wavelength/(numElements*spacing));   % first null of Sum beam
fprintf('--- SYSTEM ---\n');
fprintf('Carrier %.2f GHz | wavelength %.1f cm | spacing %.2f cm\n', ...
    carrierFreq/1e9, wavelength*100, spacing*100);
fprintf('Sum-beam first nulls: +/- %.2f deg\n\n', nullAngle);

%% 2. Signal and target ---------------------------------------------------
% TARGET: a far-field RF emitter. Far away -> its wave is a plane wave at
%         the array, described only by its direction (azimuth, elevation).
true_target_angle = 6.15;                         % azimuth (deg) - unknown to the receiver
targetEl          = 0;                            % elevation (deg), emitter in array plane
targetDir         = [true_target_angle; targetEl];

% SIGNAL: a radio pulse at the carrier (1.5 GHz). After down-conversion
%         (or direct sampling in the RFSoC) one pulse sample is a complex
%         number A*exp(j*phi). Monopulse needs only ONE pulse/snapshot.
signalAmp   = 1;                                  % amplitude (V)
signalPhase = 2*pi*rand;                          % unknown transmitter phase
s           = signalAmp*exp(1j*signalPhase);      % transmitted pulse sample
signalPower = abs(s)^2;

% NOISE: receiver thermal noise, complex Gaussian, independent on each antenna
SNRdB      = 10;                                  % per-antenna SNR
noisePower = signalPower/10^(SNRdB/10);

% RECEIVED SNAPSHOT: collectPlaneWave gives each antenna the same pulse with
% a phase shift set by the path difference d*sin(theta):  x = s*a(theta) + n
xClean = collectPlaneWave(ula, s, targetDir, carrierFreq, speed);         % 1 x 8
noise  = sqrt(noisePower/2)*(randn(1,numElements) + 1j*randn(1,numElements));
xRx    = xClean + noise;                          % what the 8 antennas deliver

fprintf('--- SIGNAL & TARGET ---\n');
fprintf('Target azimuth %.2f deg | signal power %.2f | noise power %.3f | SNR %d dB\n', ...
    true_target_angle, signalPower, noisePower, SNRdB);
fprintf('Phase step between neighbouring antennas: %.2f deg\n\n', ...
    360*spacing*sind(true_target_angle)/wavelength);

figure('Name','Received signal on each antenna','Position',[80 80 1000 450]);
subplot(1,2,1);
stem(1:numElements, abs(xClean), 'b', 'filled', 'LineWidth', 1.5); hold on; grid on;
stem(1:numElements, abs(xRx), 'r', 'LineWidth', 1.2);
xlabel('Antenna number'); ylabel('Amplitude');
legend({'Without noise','With noise'}, 'Location','south');
title('Amplitude: equal on all antennas');
subplot(1,2,2);
plot(1:numElements, unwrap(angle(xClean))*180/pi, 'b-o', 'LineWidth', 1.8); hold on; grid on;
plot(1:numElements, unwrap(angle(xRx))*180/pi,    'r--x', 'LineWidth', 1.2);
xlabel('Antenna number'); ylabel('Phase (degrees)');
legend({'Without noise','With noise'}, 'Location','northwest');
title(sprintf('Phase: linear step = direction of the %.2f^\\circ target', true_target_angle));

%% 3. Array beam: main lobe and side lobes with pattern() ----------------
% pattern() with no weights draws the array's (sum) beam at boresight.
plotAngles = -90:0.05:90;
pBeam  = pattern(ula, carrierFreq, plotAngles, 0, 'Type','powerdb');   % normalised to 0 dB
sumdB  = max(pBeam(:).', -60);

hp    = plotAngles(sumdB >= -3);   HPBW = hp(end) - hp(1);
slIdx = find(abs(plotAngles) > nullAngle + 0.5);
[PSLL, k] = max(sumdB(slIdx));     slAngle = abs(plotAngles(slIdx(k)));
fprintf('--- BEAMS (pattern) ---\n');
fprintf('HPBW %.2f deg | null-to-null %.2f deg | peak side lobe %.2f dB at +/-%.1f deg\n\n', ...
    HPBW, 2*nullAngle, PSLL, slAngle);

figure('Name','Array beam','Position',[80 80 1000 520]);
hold on; grid on; box on;
yl = [-50 4];
patch([-nullAngle nullAngle nullAngle -nullAngle], [yl(1) yl(1) yl(2) yl(2)], [1 0.93 0.75], ...
    'EdgeColor','none', 'HandleVisibility','off');
plot(plotAngles, sumdB,  'b', 'LineWidth', 2, 'HandleVisibility','off');
plot([hp(1) hp(end)], [-3 -3], 'k-', 'LineWidth', 3, 'HandleVisibility','off');
text(0, -5.5, sprintf('HPBW = %.1f^\\circ', HPBW), 'FontWeight','bold', 'HorizontalAlignment','center');
yline(PSLL, 'm--', sprintf('Peak side lobe = %.1f dB', PSLL), 'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');
xline(-nullAngle, 'k--', sprintf('First null %.1f^\\circ', -nullAngle), 'HandleVisibility','off');
xline( nullAngle, 'k--', sprintf('First null %.1f^\\circ',  nullAngle), 'HandleVisibility','off');
xline(0, 'k:', 'Boresight', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
text(0, 2.2, 'MAIN LOBE', 'HorizontalAlignment','center', 'FontWeight','bold');
text(-(nullAngle+90)/2, 2.2, 'SIDE LOBES', 'HorizontalAlignment','center', 'FontWeight','bold');
text( (nullAngle+90)/2, 2.2, 'SIDE LOBES', 'HorizontalAlignment','center', 'FontWeight','bold');
xlim([-90 90]); ylim(yl);
xlabel('Angle (degrees)'); ylabel('Normalised gain (dB)');
title('Beam of the 8x1 array: main lobe and side lobes (pattern)');

%% 4. 360 deg polar view: main, side and back lobes (pattern) -------------
% Isotropic elements radiate equally to the back, so a realistic element
% with a 20 dB front-to-back ratio (patch / dipole over ground) is used.
FBR      = 20;
azE      = -180:180;   elE = -90:90;
[AZ, EL] = meshgrid(azE, elE);
fl       = 10^(-FBR/20);
Epat     = fl + (1-fl)*((1 + cosd(EL).*cosd(AZ))/2).^2;
elem     = phased.CustomAntennaElement('AzimuthAngles', azE, 'ElevationAngles', elE, ...
               'MagnitudePattern', 20*log10(Epat), 'PhasePattern', zeros(size(Epat)));
ulaP     = phased.ULA('NumElements', numElements, 'ElementSpacing', spacing, 'Element', elem);

az360 = -180:0.25:180;
p360  = pattern(ulaP, carrierFreq, az360, 0, 'Type','powerdb');     % normalised to 0 dB
pat   = max(p360(:).', -40);
front = abs(az360) <= 90;

tmp = pat;  tmp(~front | abs(az360) < nullAngle) = -Inf;  [sl, is] = max(tmp);   % side lobe
tmp = pat;  tmp(abs(az360) < 120) = -Inf;                 [bl, ib] = max(tmp);   % back lobe

figure('Name','360 deg polar view');
pax = polaraxes;
polarplot(pax, deg2rad(az360), pat + 40, 'b', 'LineWidth', 1.6);
pax.ThetaZeroLocation = 'top';  pax.ThetaDir = 'clockwise';
pax.RLim = [0 40];  pax.RTick = 0:10:40;
pax.RTickLabel = {'-40 dB','-30 dB','-20 dB','-10 dB','0 dB'};
txt = {'FontWeight','bold'};
text(pax, 0, 43, 'main lobe', txt{:}, 'HorizontalAlignment','center');
text(pax, deg2rad(az360(is)), sl + 46, sprintf('side lobe (%.1f dB)', sl), txt{:});
text(pax, deg2rad(az360(ib)), bl + 47, sprintf('back lobes (%.1f dB)', bl), txt{:}, ...
    'HorizontalAlignment','center');
title(pax, sprintf('Beam of the 8x1 ULA (%.1f GHz): main, side and back lobes', carrierFreq/1e9));

%% 5. Monopulse behaviour: estimated vs true angle (built-in tracker) ----
% Noiseless sweep of the target over the whole field of view. Where the
% Delta/Sum ratio is outside the tracker's range it refuses (error), which
% safeTrack records as NaN (no estimate). Side-lobe targets can still give
% a ratio inside the range -> wrong (ghost) estimates.
trueAng = -90:0.1:90;
estAng  = zeros(size(trueAng));
for i = 1:numel(trueAng)
    x = receiveSignal(ula, carrierFreq, speed, trueAng(i), Inf);
    estAng(i) = safeTrack(tracker, x, steerAng);
end
% Valid range = the continuous block of estimates around boresight (0 deg).
% Estimates outside that block are GHOSTS (side-lobe targets whose phase
% difference wrapped around 2*pi and looks like a near-boresight target).
okAng = ~isnan(estAng);
[~, i0] = min(abs(trueAng));
iL = i0;  while iL > 1               && okAng(iL-1), iL = iL - 1; end
iR = i0;  while iR < numel(trueAng)  && okAng(iR+1), iR = iR + 1; end
trkRange = trueAng([iL iR]);
inRange  = false(size(trueAng));  inRange(iL:iR) = true;
isGhost  = okAng & ~inRange;

estValid = estAng;  estValid(~inRange) = NaN;
estGhost = estAng;  estGhost(~isGhost) = NaN;

fprintf('--- MONOPULSE RANGE ---\n');
fprintf('Correct estimates : %.1f to %.1f deg (main-lobe nulls +/- %.2f deg)\n', ...
    trkRange(1), trkRange(2), nullAngle);
fprintf('Max error inside that range (noiseless): %.4f deg\n', ...
    max(abs(estAng(inRange) - trueAng(inRange))));
if any(isGhost)
    gA = trueAng(isGhost & trueAng > 0);
    fprintf('GHOST estimates   : targets at +/-%.1f to +/-%.1f deg are reported as angles near boresight\n\n', ...
        min(gA), max(gA));
else
    fprintf('No ghost estimates\n\n');
end

figure('Name','Monopulse behaviour','Position',[120 60 1000 650]);
ax1 = subplot(2,1,1); hold on; grid on; box on;
patch([trkRange(1) trkRange(2) trkRange(2) trkRange(1)], [-90 -90 90 90], [0.85 1 0.85], 'EdgeColor','none');
plot(trueAng, trueAng,  'k--', 'LineWidth', 1);
plot(trueAng, estValid, 'b',   'LineWidth', 2.5);
plot(trueAng, estGhost, 'r',   'LineWidth', 2.5);
xline(-nullAngle, 'k:', 'null', 'HandleVisibility','off');  xline(nullAngle, 'k:', 'null', 'HandleVisibility','off');
text(0, 75, sprintf('VALID \\pm%.1f^\\circ', trkRange(2)), 'HorizontalAlignment','center', 'FontWeight','bold');
text(-45, 45, 'GHOSTS (side lobes)', 'HorizontalAlignment','center', 'Color','r', 'FontWeight','bold');
text( 45, 45, 'GHOSTS (side lobes)', 'HorizontalAlignment','center', 'Color','r', 'FontWeight','bold');
text(-20, -60, 'no estimate', 'HorizontalAlignment','center');
text( 20, -60, 'no estimate', 'HorizontalAlignment','center');
xlim([-90 90]); ylim([-90 90]); ylabel('Estimated angle (degrees)');
legend({'Valid range','Ideal (estimate = true)','Correct estimate','Ghost (wrong) estimate'}, ...
    'Location','southeast');
title('Built-in monopulse: correct near boresight, refuses or gives ghosts elsewhere');
ax2 = subplot(2,1,2); hold on; grid on; box on;
patch([-nullAngle nullAngle nullAngle -nullAngle], [-50 -50 2 2], [0.85 1 0.85], 'EdgeColor','none');
plot(plotAngles, sumdB, 'b', 'LineWidth', 2);
xlim([-90 90]); ylim([-50 2]);
xlabel('True target angle (degrees)'); ylabel('Beam gain (dB)');
linkaxes([ax1 ax2], 'x');

%% 6. Single target estimation --------------------------------------------
% Uses the target and received snapshot defined in section 2
est_noiseless = safeTrack(tracker, xClean, steerAng);

nPulses = 100;                                    % repeated noisy pulses
est_pulses = zeros(1, nPulses);
est_pulses(1) = safeTrack(tracker, xRx, steerAng);   % pulse from section 2
for k = 2:nPulses
    est_pulses(k) = safeTrack(tracker, ...
        receiveSignal(ula, carrierFreq, speed, true_target_angle, SNRdB), steerAng);
end

fprintf('--- TARGET (built-in tracker) ---\n');
fprintf('True target angle             : %.3f deg\n', true_target_angle);
fprintf('Noiseless estimate            : %.3f deg (error %.3f deg)\n', est_noiseless, ...
    abs(est_noiseless - true_target_angle));
fprintf('Single noisy pulse (SNR %d dB) : %.3f deg (error %.3f deg)\n', SNRdB, est_pulses(1), ...
    abs(est_pulses(1) - true_target_angle));
fprintf('Mean of %d noisy pulses       : %.3f deg (std %.3f deg, %d pulses without estimate)\n\n', ...
    nPulses, mean(est_pulses, 'omitnan'), std(est_pulses, 'omitnan'), sum(isnan(est_pulses)));

figure('Name','Single target');
histogram(est_pulses, 25); grid on; hold on;
xline(true_target_angle, 'r-', 'True angle', 'LineWidth', 2);
xlabel('Estimated angle (degrees)'); ylabel('Number of pulses');
title(sprintf('Built-in monopulse estimates of a %.2f^\\circ target (%d pulses, SNR %d dB)', ...
    true_target_angle, nPulses, SNRdB));

%% 7. Accuracy vs angle and vs SNR (Monte-Carlo) --------------------------
% RMSE is computed over the pulses that got an estimate; the percentage
% of pulses the tracker refused (out of range) is tracked separately.
nTrials  = 300;
angTest  = -30:1:30;
rmseAng  = nan(size(angTest));
noEstAng = zeros(size(angTest));
for i = 1:numel(angTest)
    e = zeros(1, nTrials);
    for k = 1:nTrials
        e(k) = safeTrack(tracker, receiveSignal(ula, carrierFreq, speed, angTest(i), SNRdB), ...
                         steerAng) - angTest(i);
    end
    rmseAng(i)  = sqrt(mean(e.^2, 'omitnan'));
    noEstAng(i) = 100*mean(isnan(e));
end

snrList = -10:5:30;
angSNR  = [0 5 9];                                % inside the tracker's +/-9.6 deg range
rmseSNR = nan(numel(angSNR), numel(snrList));
for a = 1:numel(angSNR)
    for i = 1:numel(snrList)
        e = zeros(1, nTrials);
        for k = 1:nTrials
            e(k) = safeTrack(tracker, receiveSignal(ula, carrierFreq, speed, angSNR(a), snrList(i)), ...
                             steerAng) - angSNR(a);
        end
        rmseSNR(a,i) = sqrt(mean(e.^2, 'omitnan'));
    end
end

figure('Name','Monopulse accuracy','Position',[160 60 1000 800]);
subplot(3,1,1);
semilogy(angTest, rmseAng, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 4); grid on; hold on;
xline(-nullAngle, 'k--', 'first null');  xline(nullAngle, 'k--', 'first null');
xlim([angTest(1) angTest(end)]);
ylabel('RMSE (degrees)');
title(sprintf('Built-in monopulse error vs angle (single pulse, SNR %d dB, %d trials)', SNRdB, nTrials));
subplot(3,1,2);
plot(angTest, noEstAng, 'r-o', 'LineWidth', 1.8, 'MarkerSize', 4); grid on; hold on;
xline(-nullAngle, 'k--');  xline(nullAngle, 'k--');
xlim([angTest(1) angTest(end)]); ylim([0 100]);
xlabel('True target angle (degrees)'); ylabel('No estimate (%)');
title('Pulses the tracker refused (target outside its range)');
subplot(3,1,3);
semilogy(snrList, rmseSNR, '-o', 'LineWidth', 1.8); grid on;
xlabel('Per-antenna SNR (dB)'); ylabel('RMSE (degrees)');
legend(arrayfun(@(x) sprintf('target at %d^\\circ', x), angSNR, 'UniformOutput', false));
title('Error vs SNR: best at boresight, worse towards the main-lobe edge');

fprintf('--- ACCURACY (single pulse, SNR %d dB) ---\n', SNRdB);
for a = [0 5 10 15 18 25]
    fprintf('Target at %2d deg: RMSE %.3f deg | no estimate on %5.1f%% of pulses\n', ...
        a, rmseAng(angTest == a), noEstAng(angTest == a));
end
fprintf('\n');

%% 8. Effect of frequency across 1-2 GHz ----------------------------------
fList = [1 1.5 2]*1e9;
figure('Name','Frequency effect','Position',[200 60 1000 650]);
for i = 1:numel(fList)
    trk = phased.SumDifferenceMonopulseTracker('SensorArray', ula, ...
        'PropagationSpeed', speed, 'OperatingFrequency', fList(i));
    est = zeros(size(trueAng));
    for j = 1:numel(trueAng)
        est(j) = safeTrack(trk, receiveSignal(ula, fList(i), speed, trueAng(j), Inf), steerAng);
    end
    pf = pattern(ula, fList(i), plotAngles, 0, 'Type','powerdb');
    subplot(2,1,1); hold on;
    plot(plotAngles, max(pf(:).', -50), 'LineWidth', 1.8);
    subplot(2,1,2); hold on;
    plot(trueAng, est, 'LineWidth', 1.8);
end
subplot(2,1,1); grid on; box on; xlim([-90 90]); ylim([-50 2]); ylabel('Beam gain (dB)');
legend({'1.0 GHz','1.5 GHz','2.0 GHz'}, 'Location','southeast');
title('Higher frequency -> narrower main lobe');
subplot(2,1,2); grid on; box on; xlim([-60 60]); ylim([-60 60]);
plot(trueAng, trueAng, 'k--');
xlabel('True angle (degrees)'); ylabel('Estimated angle (degrees)');
title('...so the valid monopulse range shrinks at higher frequency');

% Tracker set for 1.5 GHz but the signal arrives at 2 GHz
xMis = receiveSignal(ula, 2e9, speed, true_target_angle, Inf);
fprintf('Target %.2f deg at 2 GHz, tracker set to 1.5 GHz: %.3f deg (wrong frequency setting)\n', ...
    true_target_angle, safeTrack(tracker, xMis, steerAng));

%% =======================================================================
function est = safeTrack(tracker, x, steer)
% Built-in monopulse; returns NaN when the measured ratio is outside the
% tracker's processing range (the tracker throws an error there).
    try
        est = tracker(x, steer);
    catch ME
        if contains(ME.message, 'beyond the monopulse processing range')
            est = NaN;
        else
            rethrow(ME);
        end
    end
end

function x = receiveSignal(ula, f, c, theta, SNRdB)
% Same model as section 2: one pulse sample s = exp(j*phi) (unit power,
% random phase) from a far-field emitter at azimuth theta, received with
% collectPlaneWave, plus complex Gaussian noise (SNRdB = Inf -> no noise).
    s = exp(1j*2*pi*rand);
    x = collectPlaneWave(ula, s, [theta; 0], f, c);
    if isfinite(SNRdB)
        noisePower = 10^(-SNRdB/10);
        x = x + sqrt(noisePower/2)*(randn(size(x)) + 1j*randn(size(x)));
    end
end
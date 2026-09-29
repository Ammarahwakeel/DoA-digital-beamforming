%% monopulse_fft_figures.m
%  Monopulse direction finding using an FFT.
%
%  8 x 1 ULA, signal somewhere in 1-2 GHz, no noise.
%  The target angle is RANDOM every run - the algorithm is not told what it is.
%
%  With the beam at boresight the formula takes its simplest form:
%
%       Sigma = R + L            (R = right four elements, L = left four)
%       Delta = R - L
%       theta = asin( (2c/(pi*N*d*f)) * atan( imag(Delta/Sigma) ) )
%
%  The answer is computed that way AND with MATLAB's built-in monopulse,
%  then the two are compared.

clear; clc; close all;
rng('shuffle');                      % new random target every run

figDir = 'figs_simple';
if ~exist(figDir,'dir'), mkdir(figDir); end

%% ---------------- 1. Setup ---------------------------------------------
lightSpeed  = physconst('LightSpeed');
bandLow     = 1e9;
bandHigh    = 2e9;
sampleRate  = 5e9;
numSamples  = 1024;
numElements = 8;

spacing    = lightSpeed/(2*bandHigh);        % half wavelength at the top of the band
array      = phased.ULA('NumElements', numElements, 'ElementSpacing', spacing);
elementPos = getElementPosition(array);      % 3 x 8
y          = elementPos(2,:);                % element positions along the array

steerAngle = 0;                              % beam points straight ahead

%% ---------------- 2. A RANDOM hidden target ----------------------------
trueAngle = (rand*2 - 1)*12;                 % random angle between -12 and +12 deg
trueFreq  = bandLow + rand*(bandHigh - bandLow);   % random frequency in the band

time     = (0:numSamples-1)'/sampleRate;
delay    = -y*sind(trueAngle)/lightSpeed;
received = cos(2*pi*trueFreq*(time - delay));      % no noise

%% ---------------- 3. FFT finds the frequency and the snapshot ----------
spectrum = fft(received .* hann(numSamples), numSamples, 1);
freqAxis = (0:numSamples-1)'/numSamples*sampleRate;      % Hz

inBand = find(freqAxis >= bandLow & freqAxis <= bandHigh & freqAxis < sampleRate/2);
[~, k]    = max(sum(abs(spectrum(inBand,:)).^2, 2));
signalBin = inBand(k);
binFreq   = freqAxis(signalBin);
foundFreq = refineFreq(spectrum, signalBin, numSamples, sampleRate);

snapshot = spectrum(signalBin, :);           % 1 x 8, the array snapshot

%% ---------------- 4. Monopulse: the formula, by hand -------------------
L = sum(snapshot(1:numElements/2));          % left four elements
R = sum(snapshot(numElements/2+1:end));      % right four elements

Sigma = R + L;
Delta = R - L;
ratio = imag(Delta/Sigma);

formulaAngle = asind( (2*lightSpeed/(pi*numElements*spacing*foundFreq)) * atan(ratio) );

%% ---------------- 5. Monopulse: MATLAB's built-in ----------------------
tracker = phased.SumDifferenceMonopulseTracker('SensorArray', array, ...
    'PropagationSpeed', lightSpeed, 'OperatingFrequency', foundFreq);
builtinAngle = tryBuiltin(tracker, snapshot, steerAngle);

%% ---------------- 6. Results -------------------------------------------
beamWidth = beamwidth(array, foundFreq, 'PropagationSpeed', lightSpeed);

fprintf('\n=============== SETUP ===============\n');
fprintf('Array      : %d elements, %.1f mm apart\n', numElements, spacing*1e3);
fprintf('Beamwidth  : %.2f deg at %.3f GHz\n', beamWidth, foundFreq/1e9);
fprintf('Beam points: %d deg (no scan)\n', steerAngle);
fprintf('\n=============== HIDDEN TARGET (random) ===============\n');
fprintf('True angle     : %8.4f deg\n', trueAngle);
fprintf('True frequency : %8.4f GHz\n', trueFreq/1e9);
fprintf('\n=============== FREQUENCY FROM THE FFT ===============\n');
fprintf('Nearest bin  : %8.4f GHz   error %+8.5f GHz\n', binFreq/1e9, (binFreq-trueFreq)/1e9);
fprintf('Interpolated : %8.4f GHz   error %+8.5f GHz\n', foundFreq/1e9, (foundFreq-trueFreq)/1e9);
fprintf('\n=============== ANGLE ===============\n');
fprintf('Monopulse FORMULA  : %9.4f deg   error %+8.4f deg\n', ...
        formulaAngle, formulaAngle-trueAngle);
fprintf('Monopulse BUILT-IN : %9.4f deg   error %+8.4f deg\n', ...
        builtinAngle, builtinAngle-trueAngle);
fprintf('\nFORMULA vs BUILT-IN differ by %.6f deg\n', abs(formulaAngle-builtinAngle));
fprintf('(monopulse ratio imag(Delta/Sigma) = %+0.6f)\n', ratio);

%% ---------------- 7. Ten more random targets ---------------------------
fprintf('\n=============== 10 MORE RANDOM TARGETS ===============\n');
fprintf('  true angle     formula     built-in    difference\n');
maxDiff = 0;
for trial = 1:10
    a = (rand*2 - 1)*12;
    f = bandLow + rand*(bandHigh - bandLow);
    rx = cos(2*pi*f*(time + y*sind(a)/lightSpeed));

    sp = fft(rx .* hann(numSamples), numSamples, 1);
    [~, kk] = max(sum(abs(sp(inBand,:)).^2, 2));
    bin  = inBand(kk);
    snap = sp(bin, :);
    fFound = refineFreq(sp, bin, numSamples, sampleRate);

    Lh = sum(snap(1:numElements/2));
    Rh = sum(snap(numElements/2+1:end));
    rr = imag((Rh-Lh)/(Rh+Lh));
    fo = asind( (2*lightSpeed/(pi*numElements*spacing*fFound)) * atan(rr) );

    tk = phased.SumDifferenceMonopulseTracker('SensorArray', array, ...
        'PropagationSpeed', lightSpeed, 'OperatingFrequency', fFound);
    bi = tryBuiltin(tk, snap, steerAngle);

    if ~isnan(bi), maxDiff = max(maxDiff, abs(bi-fo)); end
    fprintf('   %8.4f   %9.4f   %9.4f    %+9.6f\n', a, fo, bi, fo-bi);
end
fprintf('\nLargest formula vs built-in difference: %.6f deg\n\n', maxDiff);

%% ---------------- 8. FIGURES -------------------------------------------
wavelength = lightSpeed/foundFreq;
plotAngles = -60:0.05:60;
steerAll   = steervec(elementPos/wavelength, plotAngles);

% ---- Figure 1: the WHOLE spectrum, with the peak marked ----------------
f1 = figure('Color','w','Position',[60 60 1250 430]);

subplot(1,3,1);
plot(time(1:120)*1e9, received(1:120,1), 'b', 'LineWidth',1.2); grid on;
xlabel('Time (ns)'); ylabel('Amplitude');
title('Antenna 1, time domain');

subplot(1,3,2);
specDb = 20*log10(abs(spectrum(:,1))/max(abs(spectrum(:,1))));
patch([bandLow bandHigh bandHigh bandLow]/1e9, [-120 -120 10 10], ...
      [0.90 0.95 1], 'EdgeColor','none'); hold on;
plot(freqAxis/1e9, specDb, 'b', 'LineWidth',1.1);
plot(foundFreq/1e9, 0, 'ro', 'MarkerSize',11, 'LineWidth',2);
xline(sampleRate/2/1e9, 'k--', 'LineWidth',1.2);
text(sampleRate/2/1e9, -110, '  Nyquist', 'FontSize',8);
grid on; ylim([-120 8]); xlim([0 sampleRate/1e9]);
xlabel('Frequency (GHz)'); ylabel('Power (dB)');
title('Whole spectrum, 0 to f_s');
legend('1-2 GHz band','spectrum','peak','Location','south');

subplot(1,3,3);
zoomBins = signalBin + (-6:6);
stem(freqAxis(zoomBins)/1e9, ...
     20*log10(abs(spectrum(zoomBins,1))/max(abs(spectrum(:,1)))), ...
     'b', 'filled', 'LineWidth',1.2); hold on;
xline(binFreq/1e9,   'g--', 'LineWidth',1.8);
xline(foundFreq/1e9, 'r-',  'LineWidth',1.8);
xline(trueFreq/1e9,  'k:',  'LineWidth',2);
grid on; xlabel('Frequency (GHz)'); ylabel('Power (dB)');
title(sprintf('Zoom: bins are %.2f MHz apart', sampleRate/numSamples/1e6));
legend('FFT bins','nearest bin','interpolated','true','Location','south');
exportgraphics(f1, fullfile(figDir,'fig1_fft_spectrum.png'), 'Resolution',150);

% ---- Figure 2: the beams and the beamwidth -----------------------------
f2 = figure('Color','w','Position',[80 80 1000 560]);
subplot(2,1,1); hold on;
beamCentres = -60:beamWidth:60;
for bc = beamCentres
    w = steervec(elementPos/wavelength, bc);
    plot(plotAngles, 20*log10(abs(w' * steerAll)/numElements), 'LineWidth',1.1);
end
grid on; ylim([-25 2]); xlim([-60 60]);
xlabel('Angle (deg)'); ylabel('Beam response (dB)');
title(sprintf('%d beams of %.2f deg beamwidth cover 120 deg', ...
      numel(beamCentres), beamWidth));

subplot(2,1,2); hold on;
w0 = steervec(elementPos/wavelength, 0);
plot(plotAngles, 20*log10(abs(w0' * steerAll)/numElements), 'b', 'LineWidth',2);
yline(-3, 'r--', 'LineWidth',1.5);
plot([-beamWidth/2 beamWidth/2], [-3 -3], 'r', 'LineWidth',3);
text(0, -1.4, sprintf('beamwidth = %.2f deg', beamWidth), ...
     'HorizontalAlignment','center','Color','r','FontWeight','bold');
grid on; ylim([-30 2]); xlim([-40 40]);
xlabel('Angle (deg)'); ylabel('Beam response (dB)');
title('One beam, with the 3 dB beamwidth marked');
exportgraphics(f2, fullfile(figDir,'fig2_beams_and_beamwidth.png'), 'Resolution',150);

% ---- Figure 3: sum and difference beams, and the answer ----------------
f3 = figure('Color','w','Position',[80 80 950 500]);
wS      = steervec(elementPos/wavelength, steerAngle);
diffTap = sign(y).';
plot(plotAngles, 20*log10(abs(wS' * steerAll)/numElements), 'b', 'LineWidth',2); hold on;
plot(plotAngles, 20*log10(abs((wS.*diffTap)' * steerAll)/numElements), 'r--', 'LineWidth',1.6);
xline(trueAngle,    'k:', 'LineWidth',2.5);
xline(formulaAngle, 'm-', 'LineWidth',1.8);
grid on; ylim([-40 2]); xlim([-40 40]);
xlabel('Angle (deg)'); ylabel('Response (dB)');
title(sprintf(['Monopulse: \\Delta nulls at the beam direction (%d deg), ' ...
               'target found at %.3f deg'], steerAngle, formulaAngle));
legend('\Sigma sum beam','\Delta difference beam','true target','monopulse answer', ...
       'Location','southwest');
exportgraphics(f3, fullfile(figDir,'fig3_sum_and_difference.png'), 'Resolution',150);

% ---- Figure 4: formula vs built-in -------------------------------------
f4 = figure('Color','w','Position',[80 80 1100 470]);
subplot(1,2,1);
zoomAngles = trueAngle + (-2:0.005:2);
plot(zoomAngles, 20*log10(abs(wS' * steervec(elementPos/wavelength, zoomAngles))/numElements), ...
     'Color',[0.6 0.6 0.6], 'LineWidth',1.5); hold on;
xline(trueAngle,    'k-',  'LineWidth',2.5);
xline(formulaAngle, 'm--', 'LineWidth',2.5);
xline(builtinAngle, 'b:',  'LineWidth',2.5);
grid on; xlim([trueAngle-2 trueAngle+2]);
xlabel('Angle (deg)'); ylabel('Sum beam (dB)');
title('Both methods land on the target');
legend('sum beam','TRUE angle','formula','built-in','Location','south');

subplot(1,2,2);
errs = max(abs([formulaAngle builtinAngle] - trueAngle), 1e-7);
b = bar(errs, 0.5, 'FaceColor','flat');
b.CData = [0.8 0.2 0.8; 0.2 0.3 0.9];
set(gca,'YScale','log');
xticklabels({'formula','built-in'});
ylabel('|error| in degrees  (log scale)'); grid on; ylim([1e-7 1]);
for q = 1:2
    text(q, errs(q)*1.8, sprintf('%.5f', errs(q)), ...
         'HorizontalAlignment','center','FontWeight','bold');
end
title('Error against the true angle');
exportgraphics(f4, fullfile(figDir,'fig4_formula_vs_builtin.png'), 'Resolution',150);

fprintf('Figures saved in %s\n\n', fullfile(pwd, figDir));

%% ---------------- local function ---------------------------------------
function f = refineFreq(spectrum, bin, numSamples, sampleRate)
% A tone almost never lands exactly on an FFT bin. Fitting a parabola
% through the peak bin and its two neighbours (in log power) finds where
% the peak really is, to a small fraction of a bin.
p  = sum(abs(spectrum).^2, 2);
yL = log(p(bin-1));  yM = log(p(bin));  yR = log(p(bin+1));
offset = 0.5*(yL - yR)/(yL - 2*yM + yR);
f = (bin - 1 + offset)/numSamples*sampleRate;
end

function a = tryBuiltin(tk, snap, steer)
% The built-in refuses targets it considers too far from the steer direction
% (roughly half a beamwidth). The hand formula has no such restriction - it
% works right up to the ambiguity limit. Return NaN when it declines.
try
    a = tk(snap, steer);
catch
    a = NaN;
end
end

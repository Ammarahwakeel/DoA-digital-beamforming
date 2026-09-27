clear; clc; close all;

%% 1. Signal & System Setup
fs = 1000;                      % Sampling frequency (Hz)
t = 0:1/fs:0.3;                 % Time vector
s = zeros(size(t));
s = s(:);
s(201:205) = s(201:205) + 1;    % A simple pulse (our test signal)

carrierFreq = 1.5e9;            % Operating frequency (1-2 GHz band)
speed = physconst('LightSpeed');
wavelength = speed/carrierFreq;

numElements = 8;
ula = phased.ULA('NumElements', numElements, 'ElementSpacing', wavelength/2);
ula.Element.FrequencyRange = [1e9 2e9];

%% 2. Simulate Received Signal (Target at +30 degrees, no noise)
inputAngle = [30; 0];
rxSignal = collectPlaneWave(ula, s, inputAngle, carrierFreq);

%% 3. Apply a Single FFT-Domain (Subband) Delay-Sum Beamformer, Steered at the Target
% Instead of applying the phase shift directly on the time samples,
% this beamformer splits the signal into frequency subbands using an
% FFT, applies the phase shift per subband, then converts back to the
% time domain (IFFT).
%
% NOTE: 'Direction' is a non-tunable property on this System object --
% once the object has been called (stepped), MATLAB locks it and won't
% let you change it without calling release() first. Since we need to
% re-steer this same beamformer to ~120 different angles in the sweep
% below, we instead set 'DirectionSource' to 'Input port'. This lets us
% pass the steering angle as an ARGUMENT on every call instead, which
% is fully tunable call-to-call and much faster than release()-ing the
% object every loop iteration.
fftBeamformer = phased.SubbandPhaseShiftBeamformer('SensorArray', ula, ...
    'OperatingFrequency', carrierFreq, ...
    'SampleRate', fs, ...
    'NumSubbands', 64, ...
    'DirectionSource', 'Input port', ...
    'WeightsOutputPort', true);

[yOut, w] = fftBeamformer(rxSignal, inputAngle);

%% 3.5 Coarse DoA Search: Sweep the Beam Across the Area of Interest (AoI)
% This is the "Phase 1" search step: instead of steering only to the
% known target angle, we steer the SAME beamformer to every candidate
% angle in the AoI, measure the output rxPower each time, and pick the
% angle that gave the strongest response as the coarse angle estimate.

theta_s_range = -60:1:60;               % Area of Interest, 1-degree steps
rxPower = zeros(size(theta_s_range));     % preallocate output rxPower for each angle

for idx_theta = 1:length(theta_s_range)
    theta_s = theta_s_range(idx_theta);

    yOut_sweep = fftBeamformer(rxSignal, [theta_s; 0]);  % steer + apply weights + sum, in one call

    rxPower(idx_theta) = sum(abs(yOut_sweep).^2);  % total output rxPower at this angle
end

[maxPower, idx] = max(rxPower);
coarse_angle = theta_s_range(idx);

fprintf('Coarse angle estimate: %d degrees (true angle was %d degrees)\n', ...
    coarse_angle, inputAngle(1));

%% 4. VISUALIZATION: Polar Beam Pattern (steered back to the true angle for the demo plot)
[yOut, w] = fftBeamformer(rxSignal, inputAngle);

%% 3.6 Inspect the Weight Vector w_n
% Because this is a SUBBAND beamformer, w is not a single 8-element
% vector -- it's an (NumElements x NumSubbands) matrix: one complete
% w_n = A_n * exp(-j*k*d_n*sin(theta_s)) vector PER frequency subband,
% since k = 2*pi*f/c is slightly different at each subband's frequency.
% size(w) reads [8 64] here: rows = elements, columns = subbands.

fprintf('\nSize of weight matrix w: [%d elements x %d subbands]\n', size(w,1), size(w,2));

% Look at the weights for one representative subband (the center one,
% closest to the carrier frequency) -- this is the cleanest single
% w_n vector to compare directly against the w_n = A_n*exp(-jkd_n*sin(theta_s))
% formula. Since columns = subbands, we pick one COLUMN, which gives us
% all numElements weights (one per antenna) for that single subband.
centerSubband = round(size(w,2)/2);
w_center = w(:, centerSubband);        % numElements x 1 column vector

fprintf('\nWeights for center subband (#%d of %d), steered to %d degrees:\n', ...
    centerSubband, size(w,2), inputAngle(1));
disp(table((1:numElements)', abs(w_center), rad2deg(angle(w_center)), ...
    'VariableNames', {'Element_n', 'Magnitude_An', 'Phase_deg'}));

fprintf('\nSame weights as raw complex numbers (real + j*imag):\n');
for n = 1:numElements
    fprintf('  w_%d = %.4f + %.4fj\n', n, real(w_center(n)), imag(w_center(n)));
end

figure;
pattern(ula, carrierFreq, -180:1:180, 0, 'PropagationSpeed', speed, ...
    'Type', 'powerdb', ...
    'CoordinateSystem', 'polar', ...
    'Weights', w);
title(sprintf('FFT-Domain Beam Pattern Steered to +%d^\\circ (1.5 GHz)', inputAngle(1)));

%% 4.5 VISUALIZATION: Power vs. Scan Angle (the sweep result)
figure;
plot(theta_s_range, rxPower, 'LineWidth', 1.5);
hold on;
plot(coarse_angle, maxPower, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Steering Angle \theta_s (degrees)');
ylabel('Output Power');
title('Beamformer Sweep: Output Power vs. Steering Angle');
legend('Power sweep', sprintf('Peak at %d^\\circ', coarse_angle));
grid on;

clear; clc; close all;

%% 1. Signal & System Setup
t = 0:0.001:0.3;                % Time, sampling frequency is 1kHz
s = zeros(size(t)); 
s = s(:);                       % Signal in column vector
s(201:205) = s(201:205) + 1;    % Define the pulse

carrierFreq = 30e6;
speed = physconst('LightSpeed');
wavelength = speed/carrierFreq;

numElements = 8;
ula = phased.ULA('NumElements', numElements, 'ElementSpacing', wavelength/2);
ula.Element.FrequencyRange = [90e5 110e6];

%% 2. Calculate the 8 Orthogonal Angles
% For an N-element ULA with lambda/2 spacing, orthogonal beams occur at:
% theta = asin(2k/N), where k ranges from -N/2 to (N/2 - 1)
k = -numElements/2 : (numElements/2 - 1);  % For N=8, k = -4 to 3
orthogonalAngles = asind((2*k)/numElements); 

%% 3. Simulate Received Signal (NO NOISE)
% Target placed at exactly +30 degrees.
inputAngle = [30; 0];

% rxSignal is now perfectly clean without the randn() noise addition
rxSignal = collectPlaneWave(ula, s, inputAngle, carrierFreq);

%% 4. Apply 8 Beamformers & Collect Outputs
allWeights = zeros(numElements, numElements);
yBeams = zeros(length(t), numElements);

for i = 1:numElements
    psbeamformer = phased.PhaseShiftBeamformer('SensorArray', ula, ...
        'OperatingFrequency', carrierFreq, ...
        'Direction', [orthogonalAngles(i); 0], ...
        'WeightsOutputPort', true);
    
    [yCbf, w] = step(psbeamformer, rxSignal);
    
    allWeights(:, i) = w;       % Store complex weights for the spatial plot
    yBeams(:, i) = yCbf;        % Store filtered time-domain signal
end

%% 5. VISUALIZATION 1: The Spatial Beams
figure(1);
pattern(ula, carrierFreq, -90:0.1:90, 0, 'PropagationSpeed', speed, ...
    'Type', 'powerdb', 'CoordinateSystem', 'rectangular', 'Weights', allWeights);
title('FYDP: 8 Orthogonal Beams for an 8x1 ULA');
xlabel('Azimuth Angle (degrees)');
ylabel('Power (dB)');
grid on; 

%% 6. VISUALIZATION 2: Time-Domain Output for Each Beam
figure(2);
for i = 1:numElements
    subplot(4, 2, i);
    plot(t, abs(yBeams(:, i))); 
    axis tight;
    ylim([0 max(abs(yBeams(:)))*1.1]); % Keep y-axis scale identical for comparison
    title(['Beam at ', num2str(round(orthogonalAngles(i),1)), '^\circ']);
    xlabel('Time (s)'); ylabel('Mag (V)');
end
sgtitle('Time-Domain Output of All 8 Beams (Clean Signal, Target at +30^\circ)');
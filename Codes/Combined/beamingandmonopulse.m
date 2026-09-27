clear; clc; close all;

%% 1. System Setup
carrierFreq = 30e6;
speed = physconst('LightSpeed');
wavelength = speed/carrierFreq;

% 8x1 ULA 
numElements = 8;
ula = phased.ULA('NumElements', numElements, 'ElementSpacing', wavelength/2);
sv = phased.SteeringVector('SensorArray', ula, 'PropagationSpeed', speed);

%% 2. Simulate the "Hidden" Target
true_target_angle = 42.84;  % mystery angle
v_mystery = sv(carrierFreq, true_target_angle); % Raw incoming wavefron

%% STEP 1: THE BEAMFORMING SEARCH (COARSE ANGLE)

scanAngles = -90:1:90; 
received_power = zeros(size(scanAngles));

for i = 1:length(scanAngles)
    % Steer the beam to the current scan angle
    w_scan = sv(carrierFreq, scanAngles(i));
    
    % Measure the signal power arriving from that direction
    received_power(i) = abs(w_scan' * v_mystery)^2;
end

% Find the single angle that returned the highest power
[~, max_idx] = max(received_power);
coarse_angle = scanAngles(max_idx); % rough estimate

fprintf('--- Step 1: Beamforming Search ---\n');
fprintf('Peak power detected at Coarse Angle: %.1f degrees\n\n', coarse_angle);

%% STEP 2: THE MONOPULSE TRACK (FINE ANGLE REFINEMENT)

lookAngle = coarse_angle;   % Lock the monopulse center to the rough estimate
squintAngle = 5;            % Squint beams 5 degrees left and right

w_left  = sv(carrierFreq, lookAngle - squintAngle);
w_right = sv(carrierFreq, lookAngle + squintAngle);

% Form Sum and Difference
w_sum = w_left + w_right;
w_delta = w_left - w_right;

% Build a highly precise "Local Lookup Table" (+/- 5 degrees around our coarse angle)
local_test_angles = (lookAngle - 5) : 0.01 : (lookAngle + 5);
v_local = sv(carrierFreq, local_test_angles);

% Generate the local Error Curve
sum_response = w_sum' * v_local;
delta_response = w_delta' * v_local;
local_monopulse_ratio = real(delta_response ./ sum_response);

% Process the actual mystery target through our focused monopulse filters
target_sum = w_sum' * v_mystery;
target_delta = w_delta' * v_mystery;
measured_ratio = real(target_delta / target_sum);

% Look up the exact angle in our local table
[~, closest_index] = min(abs(local_monopulse_ratio - measured_ratio));
exact_estimated_angle = local_test_angles(closest_index);

fprintf('--- Step 2: Monopulse Tracking ---\n');
fprintf('True Target Angle:      %.3f degrees\n', true_target_angle);
fprintf('Coarse Estimate:        %.3f degrees\n', coarse_angle);
fprintf('Exact Monopulse Angle:  %.3f degrees\n', exact_estimated_angle);
fprintf('Final Calculation Error: %.3f degrees\n', abs(true_target_angle - exact_estimated_angle));


% Figure 1: The Coarse Sweep
figure(1);
plot(scanAngles, 10*log10(received_power), 'LineWidth', 1.5);
title('Step 1: Coarse Beamforming Search');
xlabel('Scan Angle (degrees)'); ylabel('Received Power (dB)');
grid on; hold on;
xline(coarse_angle, 'r--', 'Detected Peak', 'LabelVerticalAlignment', 'bottom');
hold off;

% Figure 2: The Local Monopulse Track
figure(2);
plot(local_test_angles, local_monopulse_ratio, 'LineWidth', 2);
hold on;
plot(exact_estimated_angle, measured_ratio, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title(sprintf('Step 2: Local Monopulse Curve (Centered at %d^\\circ)', coarse_angle));
xlabel('Angle (degrees)'); ylabel('\Delta / \Sigma Ratio');
grid on; 
ylim([-1.5 1.5]); 
xline(coarse_angle, 'k--', 'Coarse Center', 'LabelVerticalAlignment', 'bottom');
text(exact_estimated_angle + 0.2, measured_ratio - 0.1, ...
    sprintf(' Exact Target: %.2f^\\circ', exact_estimated_angle), 'Color', 'r', 'FontWeight', 'bold');
hold off;

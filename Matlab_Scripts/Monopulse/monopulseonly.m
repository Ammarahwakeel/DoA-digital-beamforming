clear; clc; close all;

%% 1. System Setup (8x1 Antenna Array)
carrierFreq = 30e6;
speed = physconst('LightSpeed');
wavelength = speed/carrierFreq;

numElements = 8;
ula = phased.ULA('NumElements', numElements, 'ElementSpacing', wavelength/2);

%% 2. Generate the Error Curve (Phase-Comparison Monopulse)
% We simulate incoming signals from -20 to +20 degrees.
testAngles = -20:0.01:20;

% Get the raw incoming waveforms for all 8 antennas across all test angles
sv = phased.SteeringVector('SensorArray', ula, 'PropagationSpeed', speed);
v_all = sv(carrierFreq, testAngles); % Size: 8 antennas x 4001 angles

% =========================================================================
% STEP 2: THE SUB-ARRAY SPLIT (NO BEAM STEERING)
% We do not generate steered weights. We just physically group the hardware.
% Sub-array A: Left half (Antennas 1, 2, 3, 4)
% Sub-array B: Right half (Antennas 5, 6, 7, 8)
% =========================================================================
sub_A_responses = sum(v_all(1:4, :), 1); 
sub_B_responses = sum(v_all(5:8, :), 1);

% Pure Monopulse Math (Sum and Difference of the two halves)
sum_response = sub_A_responses + sub_B_responses;
delta_response = sub_A_responses - sub_B_responses;

% Because Phase-Comparison Monopulse relies on the physical distance 
% between the two sub-arrays, the mathematical difference is 90 degrees out 
% of phase. Therefore, we extract the imaginary part!
monopulse_ratio = imag(delta_response ./ sum_response);

%% 3. PRACTICAL EXAMPLE: TRACKING A TARGET WITHOUT BEAMFORMING
true_target_angle = 6.15; % The hidden target
v_mystery = sv(carrierFreq, true_target_angle);

% Split the 8 raw incoming antenna feeds into the two hardware sub-arrays
target_sub_A = sum(v_mystery(1:4));
target_sub_B = sum(v_mystery(5:8));

% The hardware algebraically calculates the Sum, Delta, and Ratio
target_sum = target_sub_A + target_sub_B;
target_delta = target_sub_A - target_sub_B;
measured_ratio = imag(target_delta / target_sum);

% Look up the exact angle in the Error Curve
[~, closest_idx] = min(abs(monopulse_ratio - measured_ratio));
estimated_angle = testAngles(closest_idx);

fprintf('--- STEP 2: Phase-Comparison Monopulse (8x1 Array, No Steering) ---\n');
fprintf('True Target Angle:      %.3f degrees\n', true_target_angle);
fprintf('Hardware Measured Ratio: %.4f\n', measured_ratio);
fprintf('Estimated Target Angle:  %.3f degrees\n', estimated_angle);
fprintf('Calculation Error:       %.3f degrees\n', abs(true_target_angle - estimated_angle));

%% 4. VISUALIZATION
figure(1);
plot(testAngles, monopulse_ratio, 'b', 'LineWidth', 2);
hold on;
plot(estimated_angle, measured_ratio, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('Step 2: Phase-Comparison Monopulse Error Curve (8x1 Array)');
xlabel('Actual Target Angle (degrees)');
ylabel('Monopulse Ratio (Voltage Error)');
grid on;
ylim([-3 3]); 
xline(0, 'k--', 'Boresight', 'LabelVerticalAlignment', 'bottom');
text(estimated_angle + 0.5, measured_ratio, sprintf(' Target Found at %.2f^\\circ', estimated_angle), 'Color', 'r', 'FontWeight', 'bold');
hold off;
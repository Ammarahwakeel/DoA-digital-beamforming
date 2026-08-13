clear; clc; close all;

%% 1. System Setup
carrierFreq = 30e6;
speed = physconst('LightSpeed');
wavelength = speed/carrierFreq;

% 8x1 ULA 
numElements = 8;
ula = phased.ULA('NumElements', numElements, 'ElementSpacing', wavelength/2);
sv = phased.SteeringVector('SensorArray', ula, 'PropagationSpeed', speed);

%% 2. Simulate the "Hidden" Target and ADD NOISE
% Moved the target to the edge of the radar's vision!
true_target_angle = 15.45;  
v_clean = sv(carrierFreq, true_target_angle); 

% Freeze the noise to a specific, disruptive pattern for the demo
rng(5); 
SNR_dB = 8; % 8 dB SNR creates realistic edge-case interference
noise_power = 10^(-SNR_dB/10);
noise = sqrt(noise_power/2) * (randn(numElements,1) + 1i*randn(numElements,1));

v_noisy = v_clean + noise; 

%% ========================================================================
%% TEST 1: MONOPULSE ONLY (NO BEAMFORMING)
%% ========================================================================
testAngles_unsteered = -20:0.001:20;
v_unsteered_sweep = sv(carrierFreq, testAngles_unsteered);
subA_sweep = sum(v_unsteered_sweep(1:4, :), 1);
subB_sweep = sum(v_unsteered_sweep(5:8, :), 1);
curve_mono_only = imag((subA_sweep - subB_sweep) ./ (subA_sweep + subB_sweep));

target_subA = sum(v_noisy(1:4));
target_subB = sum(v_noisy(5:8));
meas_ratio_mono = imag((target_subA - target_subB) / (target_subA + target_subB));

[~, idx_mono] = min(abs(curve_mono_only - meas_ratio_mono));
exact_angle_mono = testAngles_unsteered(idx_mono);
error_mono = abs(true_target_angle - exact_angle_mono);

%% ========================================================================
%% TEST 2: THE BEAMFORMING-INTEGRATED SEARCH & TRACK
%% ========================================================================
% --- Step 2A: Coarse Beamforming Search ---
scanAngles = -90:1:90; 
received_power = zeros(size(scanAngles));
for i = 1:length(scanAngles)
    w_scan = sv(carrierFreq, scanAngles(i));
    received_power(i) = abs(w_scan' * v_noisy)^2;
end
[~, max_idx] = max(received_power);
coarse_angle = scanAngles(max_idx); 

% --- Step 2B: Local Monopulse Track ---
lookAngle = coarse_angle;   
squintAngle = 5;            
w_left  = sv(carrierFreq, lookAngle - squintAngle);
w_right = sv(carrierFreq, lookAngle + squintAngle);
w_sum = w_left + w_right;
w_delta = w_left - w_right;

local_test_angles = (lookAngle - 5) : 0.001 : (lookAngle + 5);
v_local = sv(carrierFreq, local_test_angles);
sum_response = w_sum' * v_local;
delta_response = w_delta' * v_local;
local_monopulse_ratio = real(delta_response ./ sum_response);

target_sum = w_sum' * v_noisy;
target_delta = w_delta' * v_noisy;
measured_ratio_bf = real(target_delta / target_sum);

[~, closest_index] = min(abs(local_monopulse_ratio - measured_ratio_bf));
exact_angle_bf = local_test_angles(closest_index);
error_bf = abs(true_target_angle - exact_angle_bf);

%% ========================================================================
%% 3. RESULTS OUTPUT
%% ========================================================================
fprintf('--- NOISE IMMUNITY COMPARISON (SNR = %d dB) ---\n', SNR_dB);
fprintf('True Target Angle:      %.3f degrees\n\n', true_target_angle);

fprintf('TEST 1: Monopulse Only (No Beamforming)\n');
fprintf('Estimated Angle:        %.3f degrees\n', exact_angle_mono);
fprintf('Calculation Error:      %.3f degrees\n\n', error_mono);

fprintf('TEST 2: Beamforming-Integrated Search & Track\n');
fprintf('Coarse Search Angle:    %.3f degrees\n', coarse_angle);
fprintf('Exact Tracked Angle:    %.3f degrees\n', exact_angle_bf);
fprintf('Calculation Error:      %.3f degrees\n\n', error_bf);
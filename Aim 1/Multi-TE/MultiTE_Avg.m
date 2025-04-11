% CRLB_estimation_03A_Glu_Gln.m for estimating CRLBs for simulated spectra
% for Glutamate and Glutamine from SpinWizard
% Chathu 2025 March 10th University of Calgary
% Version 1

function [FID, SequenceDiagram] = CRLB_estimation_03A_Glu_Gln()

% Load simulation data
all_data{1,1} = load('master_data_glu.mat');
all_data{1,2} = load('master_data_gln.mat');

npoints = 2048;

%% Simulation Parameters
SimPars.maxCRLB = 40;
SimPars.sw = all_data{1,1}.master_data.handles.SpectralWidth * 1e3;
SimPars.dt = 1 / SimPars.sw;
SimPars.npzf = npoints;
SimPars.np = npoints;
SimPars.specband = [0 5];
SimPars.amp = 1;
SimPars.lw = 10;
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;
SimPars.T2 = 1 / (pi * SimPars.lw);
SimPars.R2 = -1 / SimPars.T2;
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1 / SimPars.T2_2;
SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);
noise_sw = 0.0045;
SimPars.noise = noise_sw * sqrt(SimPars.sw);
SimPars.noise = SimPars.noise / 1;
SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1.0 0.3];
SimPars.ncompounds = 2;

%% Load and generate spectra
SpecAll = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);
SpecAll_with_noise = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);

for TE1 = 1:SimPars.TE1_values
    for TE2 = 1:SimPars.TE2_values
        for comp = 1:SimPars.ncompounds
            SpecAll(:, TE1, TE2, comp) = SimPars.Conc(1,comp) .* all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:);
            SpecAll_with_noise(:, TE1, TE2, comp) = ...
                SimPars.Conc(1,comp) .* all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:) + ...
                SimPars.noise * randn(1,npoints) + 1i * SimPars.noise * randn(1,npoints);
        end
    end
end

% Time and frequency axes
t = 0:SimPars.dt:(SimPars.npzf-1)*SimPars.dt;
ff = -(SimPars.sw/2):(SimPars.sw/(SimPars.npzf-1)):(SimPars.sw/2);
freq = ff / all_data{1,1}.master_data.handles.LarmorFrequency + ...
       all_data{1,1}.master_data.handles.RFOffsetRx;

% Peak regions
delta = 0.00;
ppm_hi(1,1) = 2.45 - delta; ppm_low(1,1) = 2.25 + delta;
ppm_hi(1,2) = 2.55 - delta; ppm_low(1,2) = 2.35 + delta;

TE1_array = all_data{1,1}.master_data.TE1_array;
TE2_array = all_data{1,1}.master_data.TE2_array;

%% Select 4 (TE1, TE2) pairs where TE = 35 ms
valid_pairs = [];
for i = 1:length(TE1_array)
    for j = 1:length(TE2_array)
        if TE1_array(i) + TE2_array(j) == 35
            valid_pairs = [valid_pairs; i, j];
        end
    end
end

if size(valid_pairs, 1) < 4
    error('Not enough TE1 + TE2 = 35 ms combinations for averaging.');
end

selected_pairs = valid_pairs(1:4, :);

disp('Selected TE1 and TE2 values for averaging (in ms):');
for k = 1:size(selected_pairs,1)
    t1 = TE1_array(selected_pairs(k,1));
    t2 = TE2_array(selected_pairs(k,2));
    disp(['  Pair ' num2str(k) ': TE1 = ' num2str(t1) ' ms, TE2 = ' num2str(t2) ' ms, TE = ' num2str(t1 + t2) ' ms']);
end

% Average the FIDs for Glu and Gln
FID_avg = zeros(npoints, SimPars.ncompounds);
for idx = 1:4
    t1_idx = selected_pairs(idx,1);
    t2_idx = selected_pairs(idx,2);
    for comp = 1:SimPars.ncompounds
        FID_avg(:,comp) = FID_avg(:,comp) + SpecAll_with_noise(:,t1_idx,t2_idx,comp);
    end
end
FID_avg = FID_avg / 4;

% Use average spectra
Specs_with_noise(:,1) = FID_avg(:,1);
Specs_with_noise(:,2) = FID_avg(:,2);
Specs(:,1) = FID_avg(:,1);
Specs(:,2) = FID_avg(:,2);

% Display TE from first selected pair
TE1 = TE1_array(selected_pairs(1,1));
TE2 = TE2_array(selected_pairs(1,2));
TE = TE1 + TE2;

% Plot
figure(10); clf;
hold on;
plot(freq, real(fftshift(fft(Specs_with_noise(:,1)))));
plot(freq, real(fftshift(fft(Specs_with_noise(:,2)))));
set(gca, 'XDir', 'reverse'); xlim([0 5]);
legend('Glu', 'Gln');
title(['TE1 [5,7.5,10,12.5] TE2 [30,27.5,25,22.5] TE = 35ms']);

% Adjust font sizes for better readability
set(gca, 'FontSize', 16);           % Axis tick labels
xlabel('ppm', 'FontSize', 18);      % X-axis label
ylabel('Signal Intensity', 'FontSize', 18); % Y-axis label
title('TE1 [5,7.5,10,12.5] TE2 [30,27.5,25,22.5] TE = 35ms', 'FontSize', 18);   % Title font size
legend('Glu', 'Gln', 'FontSize', 20); % Legend font size

%% CRLB Calculation
D = zeros(SimPars.np, 2*SimPars.ncompounds);
for comp = 1:SimPars.ncompounds
    D(:,2*comp-1) = reshape(Specs(:,comp),1,[]) .* exp(SimPars.R2*t) .* exp(-SimPars.R2_2*TE*1e-3);
    D(:,2*comp) = SimPars.amp * reshape(Specs(:,comp),1,[]) .* exp(SimPars.R2*t) .* t .* exp(-SimPars.R2_2*TE*1e-3);
end

F = (1/(SimPars.noise^2)) * real(ctranspose(D)*D);
invF = inv(F);
CRLBNum = zeros(1,2*SimPars.ncompounds);
for comp = 1:2*SimPars.ncompounds
   CRLBNum(comp) = sqrt(invF(comp,comp));
end

for comp = 1:SimPars.ncompounds
   CRLBamp1D(comp,1,1) = 100 * CRLBNum(2*comp-1) / SimPars.amp;
end

disp(['Glu CRLB = ' num2str(CRLBamp1D(1,1,1)) '%']);
disp(['Gln CRLB = ' num2str(CRLBamp1D(2,1,1)) '%']);

%% SNR and Linewidth Calculation
spec_fft_glu = real(fftshift(fft(Specs_with_noise(:,1))));
spec_fft_gln = real(fftshift(fft(Specs_with_noise(:,2))));
spec_span = find(freq > SimPars.specband(1) & freq < SimPars.specband(2));
freq_trimmed = freq(spec_span);
spec_fft_glu_trimmed = spec_fft_glu(spec_span);
spec_fft_gln_trimmed = spec_fft_gln(spec_span);

signal_window_glu = (freq_trimmed > 2.25) & (freq_trimmed < 2.35);
signal_window_gln = (freq_trimmed > 2.35) & (freq_trimmed < 2.45);
noise_window = (freq_trimmed > 4.9);

max_signal_glu = max(abs(spec_fft_glu_trimmed(signal_window_glu)));
max_signal_gln = max(abs(spec_fft_gln_trimmed(signal_window_gln)));
noise_std = std(spec_fft_glu_trimmed(noise_window));

SNR_glu = max_signal_glu / noise_std;
SNR_gln = max_signal_gln / noise_std;

disp('-------------------------------------------');
disp(['SNR (Glutamate) = ' num2str(SNR_glu)]);
disp(['SNR (Glutamine) = ' num2str(SNR_gln)]);
disp('-------------------------------------------');

end
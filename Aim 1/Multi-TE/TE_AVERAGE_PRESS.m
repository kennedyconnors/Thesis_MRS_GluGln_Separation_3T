function [FID, SequenceDiagram] = TE_AVERAGE_PRESSA()

% Load FID data
all_data{1,1} = load('master_data_glu.mat');
all_data{1,2} = load('master_data_gln.mat');
npoints = 2048;

% === Simulation Parameters ===
SimPars.sw = all_data{1,1}.master_data.handles.SpectralWidth * 1e3;  % Hz
SimPars.dt = 1 / SimPars.sw;
SimPars.np = npoints;
SimPars.npzf = npoints;
SimPars.amp = 1;
SimPars.lw = 10 - all_data{1,1}.master_data.handles.SpectralLineWidth;
SimPars.T2 = 1 / (pi * SimPars.lw);
SimPars.R2 = -1 / SimPars.T2;
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1 / SimPars.T2_2;
SimPars.noise = 0.001 * sqrt(SimPars.sw);
SimPars.Conc = [1.0 0.3];  % Glu, Gln
SimPars.ncompounds = 2;

SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

% === Extract TE Arrays ===
TE1_array = all_data{1,1}.master_data.TE1_array;
TE2_array = all_data{1,1}.master_data.TE2_array;

% === Preallocate Spectra ===
SpecAll = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);
SpecAll_with_noise = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);

% === Load Data and Add Noise ===
for TE1 = 1:SimPars.TE1_values
    for TE2 = 1:SimPars.TE2_values
        for comp = 1:SimPars.ncompounds
            FID_clean = all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:)';
            SpecAll(:, TE1, TE2, comp) = SimPars.Conc(comp) * FID_clean;
            SpecAll_with_noise(:, TE1, TE2, comp) = SpecAll(:, TE1, TE2, comp) + ...
                SimPars.noise * randn(npoints,1) + 1i * SimPars.noise * randn(npoints,1);
        end
    end
end

% === Choose TE1/TE2 Indices to Average ===
TE1_indices = [6, 7, 7, 7];  
TE2_indices = [6, 6, 7, 4];  

% === Compute Average TE Values ===
avg_TE1 = mean(TE1_array(1, TE1_indices));
avg_TE2 = mean(TE2_array(1, TE2_indices));
avg_TE = avg_TE1 + avg_TE2;

% === Initialize and Average Spectra ===
Specs = zeros(npoints, SimPars.ncompounds);
Specs_with_noise = zeros(npoints, SimPars.ncompounds);

for i = 1:4
    Specs(:,1) = Specs(:,1) + SpecAll(:,TE1_indices(i),TE2_indices(i), 1);  % Glu
    Specs(:,2) = Specs(:,2) + SpecAll(:,TE1_indices(i),TE2_indices(i), 2);  % Gln
    Specs_with_noise(:,1) = Specs_with_noise(:,1) + SpecAll_with_noise(:,TE1_indices(i),TE2_indices(i), 1);
    Specs_with_noise(:,2) = Specs_with_noise(:,2) + SpecAll_with_noise(:,TE1_indices(i),TE2_indices(i), 2);
end

Specs = Specs / 4;
Specs_with_noise = Specs_with_noise / 4;

% === Time and Frequency Axes ===
t = (0:SimPars.npzf-1) * SimPars.dt;
ff = linspace(-SimPars.sw/2, SimPars.sw/2, SimPars.npzf);
freq = ff ./ all_data{1,1}.master_data.handles.LarmorFrequency + ...
           all_data{1,1}.master_data.handles.RFOffsetRx;

% === Plot Spectra ===
figure(10); clf;
plot(freq, real(fftshift(fft(Specs_with_noise(:,1)))), 'LineWidth', 1.2); hold on;
plot(freq, real(fftshift(fft(Specs_with_noise(:,2)))), 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse'); xlim([0 5]);
legend('Glu', 'Gln');
title(['Avg TE1 = ' num2str(avg_TE1, '%.1f') ' ms, Avg TE2 = ' num2str(avg_TE2, '%.1f') ...
    ' ms, Total TE = ' num2str(avg_TE, '%.1f') ' ms']);

% === CRLB Calculation ===
D = zeros(SimPars.np, 2 * SimPars.ncompounds);
for comp = 1:SimPars.ncompounds
    signal = reshape(Specs(:, comp), 1, []);
    decay = exp(SimPars.R2 * t) .* exp(-SimPars.R2_2 * avg_TE * 1e-3);
    D(:, 2*comp-1) = signal .* decay;
    D(:, 2*comp)   = signal .* decay .* t;
end

F = (1 / SimPars.noise^2) * real(D' * D);
invF = inv(F);
CRLBNum = sqrt(diag(invF))';

CRLB_Glu = 100 * CRLBNum(1) / SimPars.amp;
CRLB_Gln = 100 * CRLBNum(3) / SimPars.amp;

disp(['Glu CRLB = ' num2str(CRLB_Glu, '%.2f') ' %']);
disp(['Gln CRLB = ' num2str(CRLB_Gln, '%.2f') ' %']);

% === SNR Calculation ===
delta = 0.00;
ppm_hi = [2.45, 2.55];  % Glu, Gln
ppm_low = [2.25, 2.35];

spec_fft_glu = real(fftshift(fft(Specs_with_noise(:,1))));
spec_fft_gln = real(fftshift(fft(Specs_with_noise(:,2))));

noise_window = freq > 5;
signal_window_glu = (freq > ppm_low(1)) & (freq < ppm_hi(1));
signal_window_gln = (freq > ppm_low(2)) & (freq < ppm_hi(2));

max_signal_glu = max(spec_fft_glu(signal_window_glu));
max_signal_gln = max(spec_fft_gln(signal_window_gln));
noise_std = std(spec_fft_glu(noise_window));  % assume same noise for both

SNR_glu = max_signal_glu / noise_std;
SNR_gln = max_signal_gln / noise_std;

disp(['Glu SNR = ' num2str(SNR_glu, '%.2f')]);
disp(['Gln SNR = ' num2str(SNR_gln, '%.2f')]);

end

% CRLB_estimation_04_Glu_Gln.m for estimating CRLB's for simulated spectra
% for Glutamate and Glutamine from SpinWizard using TE-averaged PRESS
% Modified from CRLB_estimation_03C_Glu_Gln.m to implement TE-averaging
% 
% Chathu 2025 April 17th University of Calgary
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FID, SequenceDiagram] = CRLB_estimation_04_Glu_Gln()

% figure(10); clf; figure(11); clf; figure(12); clf;

% Load data
all_data{1,1} = load('master_data_glu.mat');
all_data{1,2} = load('master_data_gln.mat');
% all_data{1,3} = load('master_data_GABA.mat');
% all_data{1,4} = load('master_data_GSH.mat');
npoints = 2048;

% Simulation parameters
SimPars.specband = [0 5];           % Spectral band to consider for concatenating (ppm)

SimPars.maxCRLB = 40;               % Maximum CRLB (%) for display purposes

SimPars.sw = all_data{1,1}.master_data.handles.SpectralWidth*1e3;  % Spectral width in Hz
SimPars.dt = 1/SimPars.sw;          % Dwell time (s)
SimPars.npzf = npoints;             % Number of spectral points
SimPars.np = npoints;               % Number of spectral points

SimPars.amp = 1;                    % Relative amplitude
SimPars.lw = 10;                    % Linewidth (in Hz) must be larger than spectral linewidth ( = 4 Hz)
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;

SimPars.T2 = 1/(pi*SimPars.lw);     % T2 relaxation constant (in sec)
SimPars.R2 = -1/SimPars.T2;         % R2 relaxation rate (in Hz)

SimPars.T2_2 = 0.18;                % T2 relaxation in second dimension
SimPars.R2_2 = 1/SimPars.T2_2;      % R2 relaxation rate (in Hz)

SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

noise_sw = 0.015;
SimPars.noise = noise_sw*sqrt(SimPars.sw);  % FID noise level
SimPars.noise = SimPars.noise/1;

SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1.0 0.3];           % Metabolite concentrations
SimPars.ncompounds = 2;

% Define TE parameters for TE-averaging
% Based on Hurd et al., using 4 steps with TE increment of 10ms starting at TE=35ms
% Define TE parameters for TE-averaging
SimPars.initialTE = 35;             % Initial TE in ms
SimPars.TEincrement = 10;           % TE increment in ms
SimPars.numTEsteps = 4;             % Number of TE steps 
SimPars.TEarray = SimPars.initialTE:SimPars.TEincrement:(SimPars.initialTE+(SimPars.numTEsteps-1)*SimPars.TEincrement);
SimPars.effectiveTE = mean(SimPars.TEarray); % Effective TE is the average

% Display TE parameters
disp('--------- TE-Averaging Parameters ----------');
disp(['Initial TE = ' num2str(SimPars.initialTE) ' ms']);
disp(['TE increment = ' num2str(SimPars.TEincrement) ' ms']);
disp(['Number of steps = ' num2str(SimPars.numTEsteps)]);
disp(['TE range = ' num2str(min(SimPars.TEarray)) '-' num2str(max(SimPars.TEarray)) ' ms']);
disp(['Effective TE = ' num2str(SimPars.effectiveTE) ' ms']);
disp('-------------------------------------------');

% Initialize arrays for spectra
SpecAll = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);
SpecAll_with_noise = zeros(npoints, SimPars.TE1_values, SimPars.TE2_values, SimPars.ncompounds);

% Select a specific TE1 and TE2 index to work with (using index 6 as in original code)
TE1_indx = 6;
TE2_indx = 6;

TE1 = all_data{1,1}.master_data.TE1_array(1, TE1_indx);
TE2 = all_data{1,1}.master_data.TE2_array(1, TE2_indx);
totalTE = TE1 + TE2;

% Display the selected TE values
disp('-------------------------------------------');
disp(['Selected TE1 = ' num2str(TE1) ' ms']);
disp(['Selected TE2 = ' num2str(TE2) ' ms']);
disp(['Total selected TE = ' num2str(totalTE) ' ms']);
disp('-------------------------------------------');

% Initialize arrays to store TE-averaged data
TEaveraged_FIDs = zeros(npoints, SimPars.ncompounds);
TEaveraged_FIDs_with_noise = zeros(npoints, SimPars.ncompounds);

% Process each TE step in the TE array and average the results
for te_step = 1:SimPars.numTEsteps
    % Calculate the current TE
    currentTE = SimPars.TEarray(te_step);
    
    % For each compound, extract FID for the current TE or simulate appropriate T2 decay
    for comp = 1:SimPars.ncompounds
        % Get the FID for the base TE (TE1+TE2)
        baseFID = all_data{1,comp}.master_data.all_data_struct{TE1_indx,TE2_indx}.FID(1,:);
        
        % Apply additional T2 decay to match the current TE in the averaging scheme
        TEDifference = currentTE - totalTE;
        if TEDifference >= 0
            % Apply additional T2 decay if current TE is greater than base TE
            decayFactor = exp(SimPars.R2_2 * TEDifference * 1e-3);
           processedFID = SimPars.Conc(1,comp) * baseFID.' * decayFactor;
        else
            
        end
        
        % Add these lines before the error line (line 115)
        % Add to the average
        TEaveraged_FIDs(:, comp) = TEaveraged_FIDs(:, comp) + processedFID / SimPars.numTEsteps;
        processedFID = SimPars.Conc(1,comp) * baseFID.' * decayFactor;
        
        % Create version with noise
        noise_component = SimPars.noise*randn(1,npoints) + 1i*SimPars.noise*randn(1,npoints);
        TEaveraged_FIDs_with_noise(:, comp) = TEaveraged_FIDs_with_noise(:, comp) + (processedFID + noise_component.') / SimPars.numTEsteps;
    end
end

% Time-domain axis, dimension 1
t = 0:SimPars.dt:(SimPars.npzf-1)*SimPars.dt;

% Frequency domain axis
ff = -(SimPars.sw/2):(SimPars.sw/(SimPars.npzf-1)):(SimPars.sw/2);
freq = ff./all_data{1,1}.master_data.handles.LarmorFrequency + all_data{1,1}.master_data.handles.RFOffsetRx;
spec_span = find(freq>SimPars.specband(1,1) & freq<SimPars.specband(1,2));  % Find coordinates corresponding to frequency span of interest
new_sw = all_data{1,1}.master_data.handles.LarmorFrequency.*(SimPars.specband(1,2) - SimPars.specband(1,1));
dtnew = 1/new_sw;          % Dwell time of truncated data

% Apply T2 relaxation for CRLB calculation
T2mod = exp(SimPars.R2*t');

% Convert TE-averaged FIDs to spectra, apply T2 decay compensation, and trim to spectral band
temp1 = fftshift(fft(squeeze(TEaveraged_FIDs(:, 1))));
temp2 = fftshift(fft(squeeze(TEaveraged_FIDs(:, 2))));

temp1 = reshape(temp1(spec_span),1, []);
temp2 = reshape(temp2(spec_span),1, []);

Glu_FID_truncated = ifft(ifftshift(temp1));
Gln_FID_truncated = ifft(ifftshift(temp2));

FIDs_all_concat = reshape(Glu_FID_truncated,[],1);
FIDs_all_concat(:,2) = reshape(Gln_FID_truncated,[],1);

npoints_new = size(temp1,2);
t_new = 0:dtnew:(npoints_new-1)*dtnew;

% Calculate design matrix D for CRLB calculation
D = zeros(npoints_new,2*SimPars.ncompounds);
for comp = 1:SimPars.ncompounds
    % Derivative to amplitude
    D(:,2*comp-1) = reshape(FIDs_all_concat(:, comp),1,[]).*exp(SimPars.R2*t_new);

    % Derivative to relaxation rate
    D(:,2*comp-0) = SimPars.amp*reshape(FIDs_all_concat(:, comp),1,[]).*exp(SimPars.R2*t_new).*t_new;
end

% Calculate Fisher matrix
F = (1/(SimPars.noise*SimPars.noise))*real(ctranspose(D)*D);

% Invert Fisher matrix
invF = inv(F);

% Calculate numerical Cramer-Rao lower bounds
CRLBNum = zeros(1,2*SimPars.ncompounds);
for comp = 1:2*SimPars.ncompounds
   CRLBNum(comp) = sqrt(invF(comp,comp));
end

% Calculate percentage CRLB for amplitude
CRLBamp1D = zeros(SimPars.ncompounds, 1);
for comp = 1:SimPars.ncompounds
   CRLBamp1D(comp) = 100*CRLBNum(2*comp-1)/SimPars.amp;
end

% Display CRLB results
disp('---------- CRLB Results (TE-averaged) ----------');
disp(['Glu CRLB = ' num2str(CRLBamp1D(1)) '%']);
disp(['Gln CRLB = ' num2str(CRLBamp1D(2)) '%']);
disp('-----------------------------------------------');

%% --- SNR and Linewidth Calculations (based on frequency-domain signal) ---

% Convert TE-averaged FIDs with noise to frequency domain
spec_fft_glu = fftshift(fft(TEaveraged_FIDs_with_noise(:, 1)));
spec_fft_gln = fftshift(fft(TEaveraged_FIDs_with_noise(:, 2)));

% Trim to spectral band
spec_fft_glu_trimmed = spec_fft_glu(spec_span);
spec_fft_gln_trimmed = spec_fft_gln(spec_span);
freq_trimmed = freq(spec_span);

% Define metabolite signal regions - modified for TE-averaged spectra
% In TE-averaged spectra, Glu peak at 2.35 ppm is better resolved from Gln
signal_window_glu = (freq_trimmed > 2.30) & (freq_trimmed < 2.40);  % More specific for Glu H4
signal_window_gln = (freq_trimmed > 3.70) & (freq_trimmed < 3.80);  % Focus on Glx (includes Gln) at 3.75 ppm

% Define noise window
noise_window = (freq_trimmed > 4.9);  % Far from metabolite peaks

% Calculate peak signal amplitude
max_signal_glu = max(abs(spec_fft_glu_trimmed(signal_window_glu)));
max_signal_gln = max(abs(spec_fft_gln_trimmed(signal_window_gln)));

% Estimate noise level
noise_std = std(real(spec_fft_glu_trimmed(noise_window)));

% Compute SNRs
SNR_glu = max_signal_glu / noise_std;
SNR_gln = max_signal_gln / noise_std;

% Display SNR results
disp('----------- SNR Results (TE-averaged) -----------');
disp(['SNR (Glutamate) = ' num2str(SNR_glu, '%.2f')]);
disp(['SNR (Glutamine) = ' num2str(SNR_gln, '%.2f')]);
disp('------------------------------------------------');

% Plot the TE-averaged spectra
figure;
subplot(2,1,1);
plot(freq_trimmed, real(spec_fft_glu_trimmed));
title('TE-averaged Glutamate Spectrum');
xlabel('Frequency (ppm)');
ylabel('Amplitude');
xlim([1.5 4.0]);
set(gca, 'XDir', 'reverse');  % Reverse x-axis for standard NMR presentation

subplot(2,1,2);
plot(freq_trimmed, real(spec_fft_gln_trimmed));
title('TE-averaged Glutamine Spectrum');
xlabel('Frequency (ppm)');
ylabel('Amplitude');
xlim([1.5 4.0]);
set(gca, 'XDir', 'reverse');  % Reverse x-axis for standard NMR presentation

end

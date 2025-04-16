% in CRLB_estimation_03A_Glu_Gln, add capability to change spectral width,
% to be compabible with the hybrid version
% CRLB_estimation_03_Glu_Gln.m for estimating CRLB's for simulated spectra
% for Glutamate and Glutamine from SpinWizard

% Chathu 2025 March 10th University of Calgary
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FID, SequenceDiagram] = CRLB_estimation_03C_Glu_Gln()

% figure(10); clf; figure(11); clf; figure(12); clf;

%T2 of Glc taken to be between 100-170 from Wyss et al. 2018
% [FID, SequenceDiagram] = SpinWizard_SimulateHomonuclear_ECLIPSE_dSE(handles);

all_data{1,1} = load('master_data_glu.mat');
all_data{1,2} = load('master_data_gln.mat');
% all_data{1,3} = load('master_data_GABA.mat');
% all_data{1,4} = load('master_data_GSH.mat');
npoints = 2048;


SimPars.specband = [0 5];           %Spectral band to consider for concatenating********************
%%
SimPars.maxCRLB = 40;      % Maximum CRLB (%) for display purposes


SimPars.sw = all_data{1,1}.master_data.handles.SpectralWidth*1e3;   %spectral width in Hz
SimPars.dt = 1/SimPars.sw;          % Dwell time (s)
SimPars.npzf = npoints;                % Number of spectral points
SimPars.np = npoints;                % Number of spectral points


SimPars.amp = 1;                    % Relative amplitude
SimPars.lw = 10;                    % Linewidth (in Hz) must be larger than spectral linewidth ( = 4 Hz)
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;

SimPars.T2 = 1/(pi*SimPars.lw);     % T2 relaxation constant (in sec)
SimPars.R2 = -1/SimPars.T2;         % R2 relaxation rate (in Hz)

% SimPars.lw2 = 4;                    % Linewidth in second dimension (Hz)
% SimPars.T2_2 = 1/(pi*SimPars.lw2);  % T2 relaxation in second dimension
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1/SimPars.T2_2;      % R2 relaxation rate (in Hz)


SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

noise_sw = 0.015;
SimPars.noise = noise_sw*sqrt(SimPars.sw);  % FID noise level
SimPars.noise = SimPars.noise/1;

SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1.0 0.3];
SimPars.ncompounds = 2;

SpecAll = zeros(npoints, SimPars.TE1_values,SimPars.TE2_values,SimPars.ncompounds);
SpecAll_with_noise = zeros(npoints, SimPars.TE1_values,SimPars.TE2_values,SimPars.ncompounds);
for TE1=1:SimPars.TE1_values
    for TE2 = 1: SimPars.TE2_values
        for comp = 1:SimPars.ncompounds
            SpecAll(:, TE1, TE2, comp) = SimPars.Conc(1,comp).*all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:);
            SpecAll_with_noise(:, TE1, TE2, comp) = SimPars.Conc(1,comp).*all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:) + SimPars.noise*randn(1,npoints) + 1i*SimPars.noise*randn(1,npoints);
        end
    end
end


% Time-domain axis, dimension 1
t = 0:SimPars.dt:(SimPars.npzf-1)*SimPars.dt;


ff = -(SimPars.sw/2):(SimPars.sw/(SimPars.npzf-1)):(SimPars.sw/2);
freq = ff./all_data{1,1}.master_data.handles.LarmorFrequency + all_data{1,1}.master_data.handles.RFOffsetRx;        %***********
spec_span = find(freq>SimPars.specband(1,1) & freq<SimPars.specband(1,2));  %find coordinates on freq vector, corresponding to the frequency span of interest (in ppm) ***********
new_sw = all_data{1,1}.master_data.handles.LarmorFrequency.*(SimPars.specband(1,2) - SimPars.specband(1,1)); % ************
dtnew = 1/new_sw;          % Dwell time of truncated data



S = 0;
T2mod = exp(SimPars.R2*t');




TE1_array = all_data{1,1}.master_data.TE1_array;
TE2_array = all_data{1,1}.master_data.TE2_array;
num_TE1_values = size(TE1_array,2);



TE1_indx = 1;
TE2_indx = 1;

TE1 = TE1_array(1, 9);

TE2 = TE2_array(1, 9);
TE = TE1 + TE2;

% Display TE values
disp('-------------------------------------------');
disp(['TE1 = ' num2str(TE1) ' ms']);
disp(['TE2 = ' num2str(TE2) ' ms']);
disp(['Total TE = ' num2str(TE) ' ms']);
disp('-------------------------------------------');

FIDS_no_noise(:,1)=  SpecAll(:,TE1_indx,TE2_indx, 1);     %Glutamate spec without noise
FIDS_no_noise(:,2)=  SpecAll(:,TE1_indx,TE2_indx, 2);     %Glutamin spec without noise


temp1 = fftshift(fft(squeeze(FIDS_no_noise(:, 1).*exp(-SimPars.R2_2*TE*1e-3))));
temp2 = fftshift(fft(squeeze(FIDS_no_noise(:, 2).*exp(-SimPars.R2_2*TE*1e-3))));

temp1 = reshape(temp1(spec_span),1, []);
temp2 = reshape(temp2(spec_span),1, []);

Glu_FID_truncated = ifft(ifftshift(temp1));
Gln_FID_truncated = ifft(ifftshift(temp2));

FIDs_all_concat = reshape(Glu_FID_truncated,[],1);
FIDs_all_concat(:,2) = reshape(Gln_FID_truncated,[],1);

npoints_new = size(temp1,2);
t_new = 0:dtnew:(npoints_new-1)*dtnew;

D = zeros(npoints_new,2*SimPars.ncompounds);
for comp =1: SimPars.ncompounds

    % Derivative to amplitude
    D(:,2*comp-1) = reshape(FIDs_all_concat(:, comp),1,[]).*exp(SimPars.R2*t_new);

    % Derivative to relaxation rate
    D(:,2*comp-0) = SimPars.amp*reshape(FIDs_all_concat(:,  comp),1,[]).*exp(SimPars.R2*t_new).*t_new;
end

% Calculate Fisher matrix
F = (1/(SimPars.noise*SimPars.noise))*real(ctranspose(D)*D);

% Invert Fisher matrix
invF = inv(F);

% Calculate numerical Cramer-Rao lower bounds
CRLBNum = zeros(1,2*SimPars.ncompounds);
for comp = 1:1:2*SimPars.ncompounds;
   CRLBNum(comp) = sqrt(invF(comp,comp));
end;

for comp = 1:1:SimPars.ncompounds;
   CRLBamp1D(comp,TE1_indx,TE2_indx) = 100*CRLBNum(2*comp-1)/SimPars.amp;
end;                 
    
disp(['Glu CRLB =' num2str(CRLBamp1D(1,TE1_indx,TE2_indx)) '%']);
disp(['Gln CRLB =' num2str(CRLBamp1D(2,TE1_indx,TE2_indx)) '%']);
%%

%% --- SNR and Linewidth Calculations (based on frequency-domain signal) ---

% Extract FFT of spectra with noise for selected TE pair
spec_fft_glu = fftshift(fft(SpecAll_with_noise(:, TE1_indx, TE2_indx, 1)));
spec_fft_gln = fftshift(fft(SpecAll_with_noise(:, TE1_indx, TE2_indx, 2)));

% Define spectral span from 0–5 ppm
spec_span = find(freq > SimPars.specband(1) & freq < SimPars.specband(2));
freq_trimmed = freq(spec_span);
spec_fft_glu_trimmed = spec_fft_glu(spec_span);
spec_fft_gln_trimmed = spec_fft_gln(spec_span);

% Define metabolite signal regions based on known peak locations (Glu ≈ 2.30 ppm, Gln ≈ 2.40 ppm)
signal_window_glu = (freq_trimmed > 2.20) & (freq_trimmed < 2.45);  % Glu H4
signal_window_gln = (freq_trimmed > 2.20) & (freq_trimmed < 2.45);  % Gln H4

% Define noise window (assume noise-only region near 5 ppm)
noise_window = (freq_trimmed > 4.9);  % Ensure this is far from any known peaks

% Peak signal amplitude
max_signal_glu = max(abs(spec_fft_glu_trimmed(signal_window_glu)));
max_signal_gln = max(abs(spec_fft_gln_trimmed(signal_window_gln)));

% Estimate noise level using standard deviation in noise region
noise_std = std(real(spec_fft_glu_trimmed(noise_window)));

% Compute SNRs as peak amplitude / noise standard deviation
SNR_glu = max_signal_glu / noise_std;
SNR_gln = max_signal_gln / noise_std;

% Display results
disp('-------------------------------------------');
disp(['SNR (Glutamate) = ' num2str(SNR_glu, '%.2f')]);
disp(['SNR (Glutamine) = ' num2str(SNR_gln, '%.2f')]);
disp('-------------------------------------------');

% Plot Glu and Gln on the same axis
figure;
plot(freq_trimmed, real(spec_fft_glu_trimmed), 'LineWidth', 1.5); hold on;
plot(freq_trimmed, real(spec_fft_gln_trimmed), 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse');  % NMR convention: ppm increases to the left
xlabel('Frequency (ppm)');
ylabel('Amplitude');
title('Overlay of Glutamate and Glutamine Spectra');
legend({'Glutamate', 'Glutamine'});
xlim([1.5 4.5]);  % Adjust for visibility depending on your spectrum
grid on;
% CRLB_estimation_03_Glu_Gln.m for estimating CRLB's for simulated spectra
% for Glutamate and Glutamine from SpinWizard

% Chathu 2025 March 10th University of Calgary
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FID, SequenceDiagram] = CRLB_estimation_03A_Glu_Gln()

% figure(10); clf; figure(11); clf; figure(12); clf;

%T2 of Glc taken to be between 100-170 from Wyss et al. 2018
% [FID, SequenceDiagram] = SpinWizard_SimulateHomonuclear_ECLIPSE_dSE(handles);

all_data{1,1} = load('master_data_glu.mat');
all_data{1,2} = load('master_data_gln.mat');
% all_data{1,3} = load('master_data_GABA.mat');
% all_data{1,4} = load('master_data_GSH.mat');
npoints = 2048;

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

noise_sw = 0.0045;                          % Noise level per sqrt(spectral width)
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
freq = ff./all_data{1,1}.master_data.handles.LarmorFrequency + all_data{1,1}.master_data.handles.RFOffsetRx;

S = 0;
T2mod = exp(SimPars.R2*t');

delta = 0.00;
ppm_hi(1,1)= 2.45-delta;
ppm_low(1,1)= 2.25+delta;   %GluH4

ppm_hi(1,2)= 2.55-delta;
ppm_low(1,2)= 2.35+delta;   %GlnH4

ppm_hi(1,3)= 2.4-delta;
ppm_low(1,3)= 2.2+delta;   %GABAH4

ppm_hi(1,4)= 2.64-delta;
ppm_low(1,4)= 2.44+delta;   %GSHH4

TE1_array = all_data{1,1}.master_data.TE1_array;
TE2_array = all_data{1,1}.master_data.TE2_array;
num_TE1_values = size(TE1_array,2);

TE1_indx = 17;
TE2_indx = 17;

TE1 = TE1_array(1, TE1_indx);

TE2 = TE2_array(1, TE2_indx);
TE = TE1 + TE2;

Specs_with_noise(:,1)=  SpecAll_with_noise(:,TE1_indx,TE2_indx, 1);     %Glutamate spec with noise for visualization purposes
Specs_with_noise(:,2)=  SpecAll_with_noise(:,TE1_indx,TE2_indx, 2);     %Glutamin spec with noise for visualization purposes

Specs(:,1)=  SpecAll(:,TE1_indx,TE2_indx, 1);     %Glutamate spec without noise
Specs(:,2)=  SpecAll(:,TE1_indx,TE2_indx, 2);     %Glutamin spec without noise

figure(10); clf; 
hold on; plot(freq, real(fftshift(fft(squeeze(Specs_with_noise(:,1))))));     %plot Glutamate
hold on; plot(freq, real(fftshift(fft(squeeze(Specs_with_noise(:,2))))));     %plot Glutamate
set(gca, 'XDir', 'reverse'); xlim ([0 5]);
legend('Glu', 'Gln');
title(['TE1 = ' num2str(TE1) 'ms, TE2 = ' num2str(TE2) 'ms, TE = ' num2str(TE) 'ms.']);

D = zeros(SimPars.np,2*SimPars.ncompounds);
for comp =1: SimPars.ncompounds

    % Derivative to amplitude
    D(:,2*comp-1) = reshape(Specs(:, comp),1,[]).*exp(SimPars.R2*t).*exp(-SimPars.R2_2*TE*1e-3);

    % Derivative to relaxation rate
    D(:,2*comp-0) = SimPars.amp*reshape(Specs(:,  comp),1,[]).*exp(SimPars.R2*t).*t.*exp(-SimPars.R2_2*TE*1e-3);
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



%% --- SNR and Linewidth Calculations ---

% FFT of spectra
spec_fft_glu = real(fftshift(fft(Specs_with_noise(:,1))));
spec_fft_gln = real(fftshift(fft(Specs_with_noise(:,2))));

% --- Define signal and noise windows (adjust as needed) ---
% Noise assumed at >9.5 ppm, where no metabolite peaks are expected
noise_window = freq > 5;
signal_window_glu = (freq > ppm_low(1,1)) & (freq < ppm_hi(1,1));
signal_window_gln = (freq > ppm_low(1,2)) & (freq < ppm_hi(1,2));

% --- SNR Calculation ---
max_signal_glu = max(spec_fft_glu(signal_window_glu));
max_signal_gln = max(spec_fft_gln(signal_window_gln));
noise_std = std(spec_fft_glu(noise_window));  % Assuming same noise for both

SNR_glu = max_signal_glu / noise_std;
SNR_gln = max_signal_gln / noise_std;

disp(['Glu SNR = ' num2str(SNR_glu)]);
disp(['Gln SNR = ' num2str(SNR_gln)]);


% CRLB_estimation_05_Glu_Gln.m setup code for concatenated spectral fit, in
% 04, add monte carlo simulations

% Chathu 2025 March 10th University of Calgary
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FID, SequenceDiagram] = CRLB_estimation_05_Glu_Gln()

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

SimPars.specband = [0 5];           %Spectral band to consider for concatenating


SimPars.amp = 1;                    % Relative amplitude
SimPars.lw = 10;                    % Linewidth (in Hz) must be larger than spectral linewidth ( = 4 Hz)
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;

SimPars.T2 = 1/(pi*SimPars.lw);     % T2* relaxation constant (in sec)
SimPars.R2 = -1/SimPars.T2;         % R2* relaxation rate (in Hz)

% SimPars.lw2 = 4;                    % Linewidth in second dimension (Hz)
% SimPars.T2_2 = 1/(pi*SimPars.lw2);  % T2 relaxation in second dimension
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1/SimPars.T2_2;      % R2 relaxation rate (in Hz)


SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

% Set noise to match SNR ≈ 9.96 in single-TE case

SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1 0.3];
SimPars.ncompounds = 2;


numspec = 4;        % number of spectra to concatenate
 TE1_indx = [5,5,5,5];     % make sure number of indices in here match numspec
 TE2_indx = [1,1,1,1];     % make sure number of indices in here match numspec

%TE1_indx = [20 20 20 20];     % make sure number of indices in here match numspec
%TE2_indx = [20 20 20 20];

%SimPars.noise = SimPars.noise/sqrt(numspec);
noise_sw = 0.0045;
SimPars.noise_single = noise_sw * sqrt(SimPars.sw);    % gives ~SNR = 9.96 in single-TE
SimPars.noise = SimPars.noise_single * sqrt(numspec);  % adjust for fixed scan time (√4 = 2)

SpecAll = zeros(npoints, SimPars.TE1_values,SimPars.TE2_values,SimPars.ncompounds);
SpecAll_with_noise = zeros(npoints, SimPars.TE1_values,SimPars.TE2_values,SimPars.ncompounds);
for TE1=1:SimPars.TE1_values
    for TE2 = 1: SimPars.TE2_values
        for comp = 1:SimPars.ncompounds
            SpecAll(:, TE1, TE2, comp) = SimPars.Conc(1,comp).*all_data{1,comp}.master_data.all_data_struct{TE1,TE2}.FID(1,:);

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





TE1 = TE1_array(1, TE1_indx);

TE2 = TE2_array(1, TE2_indx);
TE = TE1 + TE2;

spec_span = find(freq>SimPars.specband(1,1) & freq<SimPars.specband(1,2));  %find coordinates on freq vector, corresponding to the frequency span of interest (in ppm)


Glu_Specs_with_noise_concat = [];         % variable that holds the spec data with noise, where data is concatenated
Gln_Specs_with_noise_concat = [];         % variable that holds the spec data with noise, where data is concatenated

Glu_Specs_no_noise_concat = [];         % variable that holds the spec data with no noise, where data is concatenated
Gln_Specs_no_noise_concat = [];         % variable that holds the spec data with no noise, where data is concatenated

Glu_Specs_for_CRLB_concat = [];           % variable that holds the spec data without noise, where data is concatenated for CRLB calculations
Gln_Specs_for_CRLB_concat = [];           % variable that holds the spec data without noise, where data is concatenated for CRLB calculations


for specIndx = 1:numspec

    
    FIDS_no_noise(:,specIndx,1)=  SpecAll(:,TE1_indx(1, specIndx),TE2_indx(1, specIndx), 1);     %Glutamate spec without noise
    FIDS_no_noise(:,specIndx,2)=  SpecAll(:,TE1_indx(1, specIndx),TE2_indx(1, specIndx), 2);     %Glutamin spec without noise
    
    temp1 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 1).*exp(SimPars.R2*t').*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3)) + (SimPars.noise*randn(npoints,1) + 1i*SimPars.noise*randn(npoints, 1))));
    temp2 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 2).*exp(SimPars.R2*t').*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3)) + (SimPars.noise*randn(npoints,1) + 1i*SimPars.noise*randn(npoints, 1))));
    Glu_Specs_with_noise_concat = [Glu_Specs_with_noise_concat, reshape(temp1(spec_span),1, [])];
    Gln_Specs_with_noise_concat = [Gln_Specs_with_noise_concat, reshape(temp2(spec_span),1, [])];
    
    temp1 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 1).*exp(SimPars.R2*t').*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3))));
    temp2 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 2).*exp(SimPars.R2*t').*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3))));
    Glu_Specs_no_noise_concat = [Glu_Specs_no_noise_concat, reshape(temp1(spec_span),1, [])];
    Gln_Specs_no_noise_concat = [Gln_Specs_no_noise_concat, reshape(temp2(spec_span),1, [])];    

    temp1 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 1).*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3))));
    temp2 = fftshift(fft(squeeze(FIDS_no_noise(:,specIndx, 2).*exp(-SimPars.R2_2*TE(1,specIndx)*1e-3))));
    Glu_Specs_for_CRLB_concat = [Glu_Specs_for_CRLB_concat, reshape(temp1(spec_span),1, [])];
    Gln_Specs_for_CRLB_concat = [Gln_Specs_for_CRLB_concat, reshape(temp2(spec_span),1, [])];

end
new_points = size(Gln_Specs_for_CRLB_concat,2);

new_freq = linspace(SimPars.specband(1,1), (SimPars.specband(1,2) - SimPars.specband(1,1))*numspec, new_points);
% new_freq = freq(spec_span);
% Specs_with_noise_conc = 

new_dt = SimPars.dt/numspec;          % new Dwell time (s)
t_new = 0:new_dt:(new_points.*new_dt - new_dt);

%convert from frequency domain to time domain (FID) for noiseless data

Glu_FID_no_noise_conc = ifft(ifftshift(Glu_Specs_for_CRLB_concat));
Gln_FID_no_noise_conc = ifft(ifftshift(Gln_Specs_for_CRLB_concat));

FIDs_all_concat = reshape(Glu_FID_no_noise_conc,[],1);
FIDs_all_concat(:,2) = reshape(Gln_FID_no_noise_conc,[],1);

FID_temp1=  SpecAll_with_noise(:,1,1, 1);     %Glutamate spec with noise for visualization purposes
spec_temp1 = real(fftshift(fft(squeeze(FID_temp1))));     %plot Glutamate

FID_temp2=  SpecAll_with_noise(:,1,1, 2);     %Glutamine spec with noise for visualization purposes
spec_temp2 = real(fftshift(fft(squeeze(FID_temp2))));     %plot Glutamine


figure(10); clf; 
hold on; plot(new_freq, real(Glu_Specs_with_noise_concat));     %plot Glutamate
hold on; plot(new_freq, real(Gln_Specs_with_noise_concat));     %plot Glutamate
set(gca, 'XDir', 'reverse');
legend('Glu', 'Gln');

figure(11); clf; 
hold on; plot(new_freq, real(Glu_Specs_no_noise_concat));     %plot Glutamate
hold on; plot(new_freq, real(Gln_Specs_no_noise_concat));     %plot Glutamate
set(gca, 'XDir', 'reverse');
legend('Glu', 'Gln');

try
    title(['TE1 = ' num2str(TE1) 'ms, TE2 = ' num2str(TE2) 'ms, TE = ' num2str(TE) 'ms.']);
catch
    ;
end
D = zeros(size(FIDs_all_concat,1),2*SimPars.ncompounds);

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
   CRLBamp1D(comp,1,1) = 100*CRLBNum(2*comp-1)/SimPars.amp;
end;                 
disp('-------------------------------------------');
disp(['TE1 = ' num2str(TE1) 'ms, TE2 = ' num2str(TE2) 'ms, TE = ' num2str(TE) 'ms.'])
disp(['Glu CRLB =' num2str(CRLBamp1D(1,1,1)) '%']);
disp(['Gln CRLB =' num2str(CRLBamp1D(2,1,1)) '%']);

%% do lsqcurve fit based monte-carlo simulation to estimate coefficient of variance

num_reps = 100;
x0 = [0.9, 0.9];     % Initial guess for parameters [amplitude1, amplitude2, offset]

MCarlo_results = zeros(num_reps, size(x0,2));

Glu_FID_no_noise_conc = ifft(ifftshift(Glu_Specs_no_noise_concat));
Gln_FID_no_noise_conc = ifft(ifftshift(Gln_Specs_no_noise_concat));
for rep_indx = 1:num_reps
    
%    noisy_spec =  fftshift(fft(SimPars.Conc(1,1).*Glu_FID_no_noise_conc + SimPars.Conc(1,2).*Gln_FID_no_noise_conc + (SimPars.noise*randn(1, new_points) + 1i*SimPars.noise*randn(1, new_points))));
   noisy_fid =  SimPars.Conc(1,1).*Glu_FID_no_noise_conc + SimPars.Conc(1,2).*Gln_FID_no_noise_conc + (SimPars.noise*randn(1, new_points) + 1i*SimPars.noise*randn(1, new_points));
   
   model_func = @(x, xdata) x(1) * real(Glu_FID_no_noise_conc) + x(2) * real(Gln_FID_no_noise_conc); % Amplitudes1*signal1 + Amplitude2*signal2
   options = optimoptions('lsqcurvefit', 'Display', 'off');
   x_opt = lsqcurvefit(model_func, x0, [], real(noisy_fid), [], [], options);
   MCarlo_results(rep_indx,:) = x_opt; 
end

Glu_mean = mean(MCarlo_results(:,1));
Gln_mean = mean(MCarlo_results(:,2));


Glu_CV = std(MCarlo_results(:,1))./Glu_mean*100;

Gln_CV = std(MCarlo_results(:,2))./Gln_mean*100;


disp('Monte Carlo simulation results:');
disp(['Glutamate concentrations = ' num2str(Glu_mean) ' +/- ' num2str(Glu_CV) '%']);
disp(['Glutamine concentrations = ' num2str(Gln_mean) ' +/- ' num2str(Gln_CV) '%']);
disp('-------------------------------------------');

% ----- SNR CALCULATION -----
% Define ppm axis for full concatenated spectrum
ppm_axis = new_freq;

% Choose Glu peak region: 2.25–2.35 ppm (you can refine this)
glu_peak_idx = find(ppm_axis > 2.25 & ppm_axis < 2.35);
noise_region_idx = find(ppm_axis > 0.5 & ppm_axis < 1.5);  % Adjust if needed

% Get real part of Glu spectrum with noise
glu_spec_real = real(Glu_Specs_with_noise_concat);

% Signal = peak height in target range
glu_peak_amplitude = max(abs(glu_spec_real(glu_peak_idx)));

% Noise = std of baseline region
noise_std = std(glu_spec_real(noise_region_idx));

% Compute SNR
SNR_Glu = glu_peak_amplitude / noise_std;

% Repeat for Gln: 2.35–2.45 ppm
gln_peak_idx = find(ppm_axis > 2.35 & ppm_axis < 2.45);
gln_spec_real = real(Gln_Specs_with_noise_concat);
gln_peak_amplitude = max(abs(gln_spec_real(gln_peak_idx)));
SNR_Gln = gln_peak_amplitude / noise_std;  % same noise baseline assumed

% Display results
disp('-------------------------------------------');
disp(['SNR (Glutamate) = ' num2str(SNR_Glu)]);
disp(['SNR (Glutamine) = ' num2str(SNR_Gln)]);
disp('-------------------------------------------');

% Labels for all evaluated TE sets
sets = {'Set A', 'Set B', 'Set C', 'Set D', 'Set E', 'Set F', 'Set G', 'Set H', 'Set I', 'Set J', 'E rvs', 'J rvs'};

% ----- CRLB Data -----
crlb_glu  = [7.24, 7.46, 8.23, 5.06, 2.29, 6.98, 7.41, 4.3, 6.57, 2.24, 2.08, 2.03];  % Glu CRLB
crlb_gln  = [28.23, 25.93, 26.76, 22.67, 8.18, 28.31, 27.86, 17.78, 25.51, 7.97, 7.32, 7.06];  % Gln CRLB

% ----- Monte Carlo CV Data -----
cv_glu    = [2.83, 11.4, 11.75, 7.1, 3.68, 11, 10.35, 6.04, 9.29, 3.35, 3.23, 2.89];  % Glu CV%
cv_gln    = [0.99, 124.93, 123.31, 89, 43.24, 123.8, 118.37, 67.89, 141.93, 36.11, 33.97, 30.32];  % Gln CV%

% ----- SNR Data -----
snr_glu   = [2.83, 3.12, 3.05, 2.23, 5.86, 2.87, 3.41, 6.14, 1.89, 5.52, 5.91, 9.211]; % Glu SNR
snr_gln   = [0.97, 1.42, 1.93, 1.46, 2.84, 2.56, 1.77, 2.49, 1.26, 2.11, 2.76, 2.08];  % Gln SNR

% ----- Plot CRLB -----
figure;
bar([crlb_glu; crlb_gln]', 'grouped');
title('CRLB Comparison for All TE Combinations', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('CRLB (%)', 'FontSize', 16);
set(gca, 'XTickLabel', sets, 'FontSize', 16,'XTickLabelRotation', 45);
legend({'Glu', 'Gln'}, 'Location', 'northeast', 'FontSize', 18);
grid on;

% ----- Plot Monte Carlo CV -----
figure;
bar([cv_glu; cv_gln]', 'grouped');
title('Monte Carlo CV (%) for All TE Combinations', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Coefficient of Variation (%)', 'FontSize', 16);
set(gca, 'XTickLabel', sets, 'FontSize', 16, 'XTickLabelRotation', 45);
legend({'Glu', 'Gln'}, 'Location', 'northeast', 'FontSize', 18);
grid on;

% ----- Plot SNR -----
figure;
bar([snr_glu; snr_gln]', 'grouped');
title('SNR Comparison for All TE Combinations', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Signal-to-Noise Ratio', 'FontSize', 16);
set(gca, 'XTickLabel', sets, 'FontSize', 16, 'XTickLabelRotation', 45);
legend({'Glu', 'Gln'}, 'Location', 'northeast', 'FontSize', 18);
grid on;
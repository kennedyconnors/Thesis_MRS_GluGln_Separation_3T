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
SimPars.lw = 15;                    % Linewidth (in Hz) must be larger than spectral linewidth ( = 4 Hz)
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;

SimPars.T2 = 1/(pi*SimPars.lw);     % T2* relaxation constant (in sec)
SimPars.R2 = -1/SimPars.T2;         % R2* relaxation rate (in Hz)

% SimPars.lw2 = 4;                    % Linewidth in second dimension (Hz)
% SimPars.T2_2 = 1/(pi*SimPars.lw2);  % T2 relaxation in second dimension
SimPars.T2_2 = 1.18;
SimPars.R2_2 = 1/SimPars.T2_2;      % R2 relaxation rate (in Hz)


SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

noise_sw = 0.008;                          % Noise level per sqrt(spectral width)
SimPars.noise = noise_sw*sqrt(SimPars.sw);  % FID noise level


SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1 1];
SimPars.ncompounds = 2;


numspec = 4;        % number of spectra to concatenate
 %TE1_indx = [35 20 15 10];     % make sure number of indices in here match numspec
 % TE2_indx = [5  20 25 30];     % make sure number of indices in here match numspec

%TE1_indx = [20 20 20 20];     % make sure number of indices in here match numspec
%TE2_indx = [20 20 20 20];

TE1_indx = [7 9 5 5];    
TE2_indx = [1 2 6 2];
SimPars.noise = SimPars.noise/sqrt(numspec);

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



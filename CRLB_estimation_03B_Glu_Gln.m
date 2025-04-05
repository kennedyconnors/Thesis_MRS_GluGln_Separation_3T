% CRLB_estimation_03B_Glu_Gln.m for estimating CRLB's for simulated spectra
% for Glutamate and Glutamine from SpinWizard, plot for an array of values

% Chathu 2025 March 10th University of Calgary
% Version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FID, SequenceDiagram] = CRLB_estimation_03B_Glu_Gln()
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
SimPars.lw = 8;                    % Linewidth (in Hz) must be larger than spectral linewidth ( = 4 Hz)
SimPars.lw = SimPars.lw - all_data{1,1}.master_data.handles.SpectralLineWidth;

SimPars.T2 = 1/(pi*SimPars.lw);     % T2 relaxation constant (in sec)
SimPars.R2 = -1/SimPars.T2;         % R2 relaxation rate (in Hz)

% SimPars.lw2 = 4;                    % Linewidth in second dimension (Hz)
% SimPars.T2_2 = 1/(pi*SimPars.lw2);  % T2 relaxation in second dimension
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1/SimPars.T2_2;      % R2 relaxation rate (in Hz)


SimPars.TE1_values = size(all_data{1,1}.master_data.TE1_array,2);
SimPars.TE2_values = size(all_data{1,1}.master_data.TE2_array,2);

noise_sw = 0.0045 *sqrt(150);                          % Noise level per sqrt(spectral width)
SimPars.noise = noise_sw*sqrt(SimPars.sw);  % FID noise level
SimPars.noise = SimPars.noise/1;

SimPars.ConcName = {'Glu','Gln'};
SimPars.Conc = [1.0 1.0];
SimPars.ncompounds = 2;

SpecAll = zeros(npoints, SimPars.TE1_values,SimPars.TE2_values,SimPars.ncompounds);

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

% ppm_hi= 2.4;
% ppm_low= 2.18;   %GABH4

% linewidth = 6;      % extra linebroadning. 


TE1_array = all_data{1,1}.master_data.TE1_array;
TE2_array = all_data{1,1}.master_data.TE2_array;
num_TE1_values = size(TE1_array,2);


T2 = 210;       % in ms assume all metabolites have this T2 relaxation


for TE1_indx = 1:num_TE1_values
    TE1 = TE1_array(1, TE1_indx);
    for TE2_indx =1 :num_TE1_values
        TE2 = TE2_array(1, TE2_indx);
        TE = TE1 + TE2;
        D = zeros(SimPars.np,2*SimPars.ncompounds);
        for comp =1: SimPars.ncompounds

            % Derivative to amplitude
            D(:,2*comp-1) = reshape(SpecAll(:, TE1_indx, TE2_indx, comp),1,[]).*exp(SimPars.R2*t).*exp(-SimPars.R2_2*TE*1e-3);

            % Derivative to relaxation rate
            D(:,2*comp-0) = SimPars.amp*reshape(SpecAll(:, TE1_indx, TE2_indx, comp),1,[]).*exp(SimPars.R2*t).*t.*exp(-SimPars.R2_2*TE*1e-3);
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

    end       
end

%
figure(11); 

subplot(1,2,1); 
imshow(log10(squeeze(CRLBamp1D(1,:,:))), [1.5 2.8], ...
    'YData', [TE1_array(1) TE1_array(end)], ...
    'XData', [TE2_array(1) TE2_array(end)], ...
    'InitialMagnification', 'fit'); 
colormap('jet'); 
set(gca,'YDir','normal');
axis on; 
title('log_{10}(Glu CRLB)');

subplot(1,2,2); 
imshow(log10(squeeze(CRLBamp1D(2,:,:))), [1.5 2.8], ...
    'YData', [TE1_array(1) TE1_array(end)], ...
    'XData', [TE2_array(1) TE2_array(end)], ...
    'InitialMagnification', 'fit'); 
colormap('jet'); 
set(gca,'YDir','normal');
axis on; 
title('log_{10}(Gln CRLB)');

%%%%

% Exhaustive search to determine Echo time averaging 
% Get all possible TE1-TE2 index pairs
TE1_inds = 1:SimPars.TE1_values;
TE2_inds = 1:SimPars.TE2_values;

all_pairs = [];
for i = 1:length(TE1_inds)
    for j = 1:length(TE2_inds)
        all_pairs = [all_pairs; TE1_inds(i), TE2_inds(j)];
    end
end

nPairs = size(all_pairs,1);
% Number of random 4-TE combinations to test
nSamples = 50000;  % You can increase or decrease this

% Random sampling of unique 4-combinations
combs = zeros(nSamples, 4);
for i = 1:nSamples
    combs(i,:) = randperm(nPairs, 4);  % Randomly pick 4 unique indices
end
nCombs = size(combs,1);

bestCRLB_Glu = Inf;
bestCRLB_Gln = Inf;
bestCRLB_Avg = Inf;

% To store best combinations
bestCombo_Glu = [];
bestCombo_Gln = [];
bestCombo_Avg = [];

fprintf('Testing %d combinations of 4 echo times...\n', nCombs);

for idx = 1:nCombs
    inds = combs(idx,:);
    combo = all_pairs(inds, :);  % 4 x 2 array of TE1/TE2 indices

     % Check for TE < 20ms in any of the 4 combinations
    skipCombo = false;
    for c = 1:4
        TE1_idx = combo(c,1);
        TE2_idx = combo(c,2);
        TE = TE1_array(1,TE1_idx) + TE2_array(1,TE2_idx);
        if TE < 20
            skipCombo = true;
            break;
        end
    end
    if skipCombo
        continue;  % Skip this 4-combo if any TE < 20 ms
    end
    % Initialize averaged signal matrix
    AvgSpec = zeros(SimPars.np, SimPars.ncompounds);


    for c = 1:4
        TE1_idx = combo(c,1);
        TE2_idx = combo(c,2);
        TE = TE1_array(1,TE1_idx) + TE2_array(1,TE2_idx);

        for comp = 1:SimPars.ncompounds
            signal = reshape(SpecAll(:, TE1_idx, TE2_idx, comp),1,[]);
            decay = exp(SimPars.R2 * t) .* exp(-SimPars.R2_2 * TE * 1e-3);
            AvgSpec(:, comp) = AvgSpec(:, comp) + (signal .* decay)';
        end
    end

    AvgSpec = AvgSpec / 4;  % Average over 4 combinations

    % Build derivative matrix D
    D = zeros(SimPars.np, 2*SimPars.ncompounds);
    for comp = 1:SimPars.ncompounds
        D(:,2*comp-1) = AvgSpec(:,comp);                     % Derivative wrt amplitude
        D(:,2*comp)   = AvgSpec(:,comp) .* t';               % Derivative wrt relaxation
    end

    % Fisher matrix
    F = (1/(SimPars.noise^2)) * real(ctranspose(D) * D);
    if rcond(F) < 1e-12
        continue  % Skip ill-conditioned matrices
    end

    invF = inv(F);
    CRLBNum = sqrt(diag(invF))';

    CRLB_Glu = 100 * CRLBNum(1) / SimPars.amp;
    CRLB_Gln = 100 * CRLBNum(3) / SimPars.amp;
    CRLB_avg = (CRLB_Glu + CRLB_Gln) / 2;

    if CRLB_Glu < bestCRLB_Glu
        bestCRLB_Glu = CRLB_Glu;
        bestCombo_Glu = combo;
    end

    if CRLB_Gln < bestCRLB_Gln
        bestCRLB_Gln = CRLB_Gln;
        bestCombo_Gln = combo;
    end

    if CRLB_avg < bestCRLB_Avg
        bestCRLB_Avg = CRLB_avg;
        bestCombo_Avg = combo;
    end
end

% Print best combinations
fprintf('\nBest combination for Glu CRLB: %.2f%%\n', bestCRLB_Glu);
disp(bestCombo_Glu);

fprintf('\nBest combination for Gln CRLB: %.2f%%\n', bestCRLB_Gln);
disp(bestCombo_Gln);

fprintf('\nBest combination for average CRLB: %.2f%%\n', bestCRLB_Avg);
disp(bestCombo_Avg);
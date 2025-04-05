function CRLB_4Combo_Averaging_GluGln()
  
%% === Load Data ===
glu_data = load('master_data_glu.mat');
gln_data = load('master_data_gln.mat');

TE1_array = glu_data.master_data.TE1_array;
TE2_array = glu_data.master_data.TE2_array;
npoints = 2048;

%% === Simulation Parameters ===
SimPars.sw = glu_data.master_data.handles.SpectralWidth * 1e3;
SimPars.dt = 1 / SimPars.sw;
SimPars.np = npoints;
SimPars.amp = 1;
SimPars.lw = 10 - glu_data.master_data.handles.SpectralLineWidth;
SimPars.T2 = 1 / (pi * SimPars.lw);
SimPars.R2 = -1 / SimPars.T2;
SimPars.T2_2 = 0.18;
SimPars.R2_2 = 1 / SimPars.T2_2;
SimPars.noise = 0.0045 * sqrt(SimPars.sw);
SimPars.Conc = [1.0 0.3];  % Glu, Gln
SimPars.ncompounds = 2;

t = (0:SimPars.np-1) * SimPars.dt;

%% === Precompute Spectra ===
n_TE1 = size(TE1_array, 2);
n_TE2 = size(TE2_array, 2);

SpecAll = zeros(npoints, n_TE1, n_TE2, SimPars.ncompounds);
for TE1 = 1:n_TE1
    for TE2 = 1:n_TE2
        SpecAll(:, TE1, TE2, 1) = SimPars.Conc(1) * glu_data.master_data.all_data_struct{TE1, TE2}.FID(1,:);
        SpecAll(:, TE1, TE2, 2) = SimPars.Conc(2) * gln_data.master_data.all_data_struct{TE1, TE2}.FID(1,:);
    end
end

%% === Identify Valid TE1/TE2 Index Pairs (TE ≥ 20 ms) ===
min_TE = 20;  % ms
valid_pairs = [];
for i = 1:n_TE1
    for j = 1:n_TE2
        if TE1_array(1,i) + TE2_array(1,j) >= min_TE
            valid_pairs = [valid_pairs; i, j];
        end
    end
end

nPairs = size(valid_pairs, 1);
fprintf('\n🔍 Found %d valid TE1/TE2 combinations with TE ≥ %d ms.\n', nPairs, min_TE);

%% === LIMIT: Truncate if too many combinations ===
maxPairs = 50;  % Change this to 100, 75, etc. if needed
if nPairs > maxPairs
    warning('Too many valid pairs — limiting to first %d for exhaustive search.', maxPairs);
    valid_pairs = valid_pairs(1:maxPairs, :);
end

%% === Generate All 4-Combinations from Valid Pairs ===
combo_indices = nchoosek(1:size(valid_pairs,1), 4);
nCombos = size(combo_indices, 1);
fprintf('🔬 Exhaustively evaluating %d 4-combinations...\n', nCombos);

%% === Optimization Loop ===
bestCRLB_avg = Inf;

for idx = 1:nCombos
    combo = valid_pairs(combo_indices(idx,:), :);  % 4 x 2

    AvgSpec = zeros(npoints, SimPars.ncompounds);
    valid = true;

    for k = 1:4
        i = combo(k,1);
        j = combo(k,2);
        TE1 = TE1_array(1,i);
        TE2 = TE2_array(1,j);
        TE = TE1 + TE2;

        if TE < min_TE
            valid = false;
            break;
        end

        for comp = 1:SimPars.ncompounds
            signal = reshape(SpecAll(:, i, j, comp), 1, []);
            decay = exp(SimPars.R2 * t) .* exp(-SimPars.R2_2 * TE * 1e-3);
            AvgSpec(:, comp) = AvgSpec(:, comp) + (signal .* decay)';
        end
    end

    if ~valid
        continue;
    end

    AvgSpec = AvgSpec / 4;

    %% === CRLB Estimation ===
    D = zeros(npoints, 2 * SimPars.ncompounds);
    for comp = 1:SimPars.ncompounds
        D(:,2*comp-1) = AvgSpec(:,comp);           % dS/dAmplitude
        D(:,2*comp)   = AvgSpec(:,comp) .* t';     % dS/dRelaxation
    end

    F = (1 / (SimPars.noise^2)) * real(D' * D);
    if rcond(F) < 1e-12
        continue;  % Skip ill-conditioned matrix
    end

    invF = inv(F);
    CRLBNum = sqrt(diag(invF))';
    CRLB_Glu = 100 * CRLBNum(1) / SimPars.amp;
    CRLB_Gln = 100 * CRLBNum(3) / SimPars.amp;
    CRLB_avg = (CRLB_Glu + CRLB_Gln) / 2;

    if CRLB_avg < bestCRLB_avg
        bestCRLB_Glu = CRLB_Glu;
        bestCRLB_Gln = CRLB_Gln;
        bestCRLB_avg = CRLB_avg;
        bestCombo = combo;
    end
end

%% === Output Best Result ===
fprintf('\n✅ Best 4-TE combination (TE1 + TE2 ≥ %d ms):\n', min_TE);
disp('   TE1_idx   TE2_idx   TE1_val   TE2_val   Total TE');
for k = 1:4
    i = bestCombo(k,1);
    j = bestCombo(k,2);
    te1_val = TE1_array(1,i);
    te2_val = TE2_array(1,j);
    fprintf('%10d %10d %10.2f %10.2f %10.2f\n', i, j, te1_val, te2_val, te1_val + te2_val);
end

fprintf('\nBest Glu CRLB: %.2f%%\n', bestCRLB_Glu);
fprintf('Best Gln CRLB: %.2f%%\n', bestCRLB_Gln);
fprintf('Average CRLB : %.2f%%\n', bestCRLB_avg);

end
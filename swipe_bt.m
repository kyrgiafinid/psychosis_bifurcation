clear;
clc;
close all;

%% ============================================================
%  BT SWEEP DRIVER
%  Runs network_simulation.m for a range of bt values and
%  collects the resulting peak LFP frequency for each.
%  ============================================================

% Edit this list to whatever bt values you want to scan.
bt_values = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40 50 60 70 80 90];

nB = numel(bt_values);
peakFreq    = zeros(nB,1);
peakPower   = zeros(nB,1);
rate_PY     = zeros(nB,1);
rate_IN     = zeros(nB,1);
meanGABA    = zeros(nB,1);
meanGABA_PY = zeros(nB,1);
meanGABA_IN = zeros(nB,1);

doPlots = false;   % set true to also get the full diagnostic figures per run

h = waitbar(0, 'Running simulations...');
for i = 1:nB
    fprintf('\n===== bt = %g  (%d/%d) =====\n', bt_values(i), i, nB);
    [peakFreq(i), peakPower(i), rate_PY(i), rate_IN(i), ...
        meanGABA(i), meanGABA_PY(i), meanGABA_IN(i)] = ...
        network_simulation(bt_values(i), doPlots);
    waitbar(i/nB, h, sprintf('Iteration %d / %d', i, nB));
end

%% ============================================================
%  RESULTS
%  ============================================================

results = table(bt_values(:), peakFreq, peakPower, rate_PY, rate_IN, ...
    meanGABA, meanGABA_PY, meanGABA_IN, ...
    'VariableNames', {'bt','peakFreq_Hz','peakPower','rate_PY_Hz','rate_IN_Hz', ...
    'meanGABA','meanGABA_PY','meanGABA_IN'});
disp(results);

writetable(results, 'bt_sweep_results.csv');
fprintf('\nResults saved to bt_sweep_results.csv\n');

figure('Name','bt sweep: peak frequency');
plot(bt_values, peakFreq, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('bt');
ylabel('Peak LFP frequency [Hz]');
title('Peak frequency vs bt');
grid on;

figure('Name','bt sweep: GABA conductance');
plot(bt_values, meanGABA,    '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'All neurons'); hold on;
plot(bt_values, meanGABA_PY, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'PY (onto PY)');
plot(bt_values, meanGABA_IN, '-^', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'IN (onto IN)');
hold off;
xlabel('bt');
ylabel('Mean gGABA (steady state)');
title('Steady-state GABA conductance vs bt');
legend('Location','best');
grid on;

%% ============================================================
%  GABA CONDUCTANCE vs PEAK (GAMMA) FREQUENCY
%  ============================================================

figure('Name','GABA vs peak frequency (colored by bt)');
scatter(meanGABA, peakFreq, 60, bt_values, 'filled');
c = colorbar;
c.Label.String = 'bt';
colormap(parula);
xlabel('Mean gGABA (steady state, all neurons)');
ylabel('Peak LFP frequency [Hz]');
title('GABA conductance vs peak (gamma) frequency');
grid on;
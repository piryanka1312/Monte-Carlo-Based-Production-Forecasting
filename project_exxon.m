%% =============================================================
%  PRODUCTION FORECAST — MONTE CARLO NODAL ANALYSIS (Well A)
%  Depletion-drive pool, water production = 0
% =============================================================
%
%  OVERVIEW:
%  1) Field-level material balance context (pressure/GOR/Np trends)
%  2) PI trend line (J/Ji vs pressure) fitted from multi-well data
%  3) Monte Carlo sampling of reservoir pressure & PI (4 depletion stages)
%  4) Monte Carlo VLP curve generation (test-data uncertainty)
%  5) IPR-VLP intersection (nodal analysis) per simulation
%  6) P10/P50/P90 statistics and plots
%
% =============================================================

clear; clc; close all;


OIIP_entire_pool = 10e6;   % bbl, stock-tank oil initially in place
pavg_current      = 1410;  % psig, current average reservoir pressure
oil_produced_current = 669000; % bbl, cumulative oil produced so far

oil   = [0, 0.031, 0.353, 0.669];        % cumulative oil, MMbbl
p_avg = [2350, 2100, 1720, 1410];        % corresponding avg pressure, psig

Np  = [0.120, 0.229, 0.307, 0.402, 0.471, 0.533, 0.565, 0.602, 0.641, 0.669];
GOR = [209, 208, 214, 220, 242, 240, 255, 298, 353, 365];

figure
yyaxis left
plot(oil, p_avg, '-o', 'LineWidth', 2)
ylabel('Average Pressure (psig)')

yyaxis right
plot(Np, GOR, '-s', 'LineWidth', 2)
ylabel('GOR (scf/stb)')
hold on

gorFitDegree = 3;
gorFitCoef   = polyfit(Np, GOR, gorFitDegree);
Np_smooth    = linspace(min(Np), max(Np), 200);
GOR_fit      = polyval(gorFitCoef, Np_smooth);
plot(Np_smooth, GOR_fit, 'r--', 'LineWidth', 2)

xlabel('Cumulative Oil Production (MMbbl)')
title('Pressure and GOR vs Cumulative Oil')
grid on
legend('Pressure', 'GOR', 'GOR Fit', 'Location', 'best')


Px_all = [2350, 1820, 1710, 1420, 2100, 1730, 1550, 2100, 1660, 1400, 2100, 1770, 1420];
jy_all = [1, 0.86, 0.41, 0.64, 1, 1.17, 0.83, 1, 0.79, 0.63, 1, 0.82, 0.64];

piTrendCoef = polyfit(Px_all, jy_all, 1);
Px_fit = linspace(min(Px_all), max(Px_all), 200);
jy_fit = polyval(piTrendCoef, Px_fit);

figure
hold on
plot(Px_fit, jy_fit, 'k-', 'LineWidth', 2, 'DisplayName', 'Single Trend Line')
plot(Px_all, jy_all, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Data Points')
set(gca, 'XDir', 'reverse')
xlabel('Average static pressure at datum, psig')
ylabel('J / J_i')
title('Production Forecast: Single PI Trend')
grid on
legend('Location', 'best')
ylim([0 1.5])
hold off

 
pressure = [1500, 1410, 1300, 1200, 1000, 800, 600, 400, 300, 200];
jy_pressure = polyval(piTrendCoef, pressure);

disp(table(pressure', jy_pressure', 'VariableNames', {'Pressure_psig', 'J_over_Ji'}))

j  = jy_pressure;
JA = j .* 0.22;   % Well A base PI trend
JB = j .* 0.06;   % Well B 
JC = j .* 0.19;   % Well C 
JD = j .* 0.11;   % Well D 

%  Monte Carlo setup

rng(123);            % reproducibility
Nsim = 5000;          % number of Monte Carlo simulations

Pressure_std = 20;    % psi, reservoir pressure measurement uncertainty


PI_residual_std = 0.05;  
GOR_std = 0.05;

Pressure_MC = zeros(Nsim, length(pressure));
JA_MC = zeros(Nsim, length(JA));
JB_MC = zeros(Nsim, length(JB));
JC_MC = zeros(Nsim, length(JC));
JD_MC = zeros(Nsim, length(JD));
max_OilRate_MC = zeros(Nsim, length(pressure));

for sim = 1:Nsim
    
    Pressure_MC(sim, :) = pressure + Pressure_std * randn(size(pressure));
    Pressure_MC(sim, Pressure_MC(sim, :) < 0) = 0;

    
    j_sim = polyval(piTrendCoef, Pressure_MC(sim, :));

    JA_MC(sim, :) = (j_sim .* 0.22) .* (1 + PI_residual_std * randn(size(JA)));
    JB_MC(sim, :) = (j_sim .* 0.06) .* (1 + PI_residual_std * randn(size(JB)));
    JC_MC(sim, :) = (j_sim .* 0.19) .* (1 + PI_residual_std * randn(size(JC)));
    JD_MC(sim, :) = (j_sim .* 0.11) .* (1 + PI_residual_std * randn(size(JD)));

    JA_MC(sim, JA_MC(sim, :) < 0) = 0;
    JB_MC(sim, JB_MC(sim, :) < 0) = 0;
    JC_MC(sim, JC_MC(sim, :) < 0) = 0;
    JD_MC(sim, JD_MC(sim, :) < 0) = 0;

    % Simple linear rate estimate (Q = J * P), used only as a first-pass
    % sanity-check distribution -- NOT the nodal analysis result.
    max_OilRate_MC(sim, :) = JA_MC(sim, :) .* Pressure_MC(sim, :);
end


figure
hold on
for sim = 1:100
    plot(Pressure_MC(sim, :), JA_MC(sim, :))
end
set(gca, 'XDir', 'reverse')
xlabel('Pressure (psig, resampled per simulation)')
ylabel('PI')
title('Monte Carlo Productivity Index Curves (Trend-Coupled to Pressure)')
grid on

figure
histogram(max_OilRate_MC(:, 1), 40)
xlabel('Oil Rate (linear J*P estimate)')
ylabel('Frequency')
title('Sanity-Check Distribution of Initial Oil Rate (Linear PI Model)')
grid on

P10_linear = prctile(max_OilRate_MC(:, 1), 90);
P50_linear = prctile(max_OilRate_MC(:, 1), 50);
P90_linear = prctile(max_OilRate_MC(:, 1), 10);
fprintf('\n===== Linear PI*P Sanity Check (Stage 1) =====\n');
fprintf('P10 = %.2f | P50 = %.2f | P90 = %.2f\n', P10_linear, P50_linear, P90_linear);


%  Nodal analysis setup (Vogel IPR + multi-rate VLP test data)

pressure2 = [1500 1410 1300 1200 1000 800 600 400 300 200 0];  % Pwf sweep for plotting

Q_test = [50 100 200];          % multi-rate test rates
pwf1 = [990 850 810];           % Pwf at Q_test, Stage 1 (Pr ~ 1500 psig)
pwf2 = [910 790 710];           % Stage 2 (Pr ~ 1410 psig)
pwf3 = [830 730 650];           % Stage 3 (Pr ~ 1300 psig)
pwf4 = [780 690 610];           % Stage 4 (Pr ~ 1200 psig)

VLP_std = 15;   % psi, uncertainty in multi-rate test pressure readings

nStages = 4;
OperatingRate     = zeros(Nsim, nStages);
OperatingPressure = zeros(Nsim, nStages);

pwf_baseline = {pwf1, pwf2, pwf3, pwf4};


% Monte Carlo IPR-VLP intersection (the actual solve)

for sim = 1:Nsim

    poly = cell(1, nStages);
    for stage = 1:nStages
        pwf_rand = pwf_baseline{stage} + randn(size(pwf_baseline{stage})) * VLP_std;
        poly{stage} = polyfit(Q_test, pwf_rand, 2);
    end

    for stage = 1:nStages
        Pr_k   = Pressure_MC(sim, stage);
        Qmax_k = JA_MC(sim, stage) * Pr_k / 1.8;   

        coef = poly{stage};

        Q_grid    = linspace(0.001, Qmax_k, 500);
        Pwf_vogel = Pr_k * 0.125 * (-1 + sqrt(81 - 80 * (Q_grid / Qmax_k)));
        Pwf_vlp   = polyval(coef, Q_grid);

        diff_curve = Pwf_vogel - Pwf_vlp;
        idx = find(diff(sign(diff_curve)) ~= 0, 1);   % first crossing only

        if isempty(idx)
            continue
        end

        Q1 = Q_grid(idx);     Q2 = Q_grid(idx + 1);
        d1 = diff_curve(idx); d2 = diff_curve(idx + 1);
        q    = Q1 - d1 * (Q2 - Q1) / (d2 - d1);
        p_op = polyval(coef, q);

        OperatingRate(sim, stage)     = q;
        OperatingPressure(sim, stage) = p_op;
    end
end


% Statistics and plots

PressureCase = ["Pr = 1500 psig"; "Pr = 1410 psig"; "Pr = 1300 psig"; "Pr = 1200 psig"];

MeanRate = zeros(1, nStages);
StdRate  = zeros(1, nStages);
P10 = zeros(1, nStages);
P50 = zeros(1, nStages);
P90 = zeros(1, nStages);
CI_lower = zeros(1, nStages);
CI_upper = zeros(1, nStages);
ProbabilityAbove100 = zeros(1, nStages);
ProbabilityAbove80  = zeros(1, nStages);

for stage = 1:nStages
    Rate = OperatingRate(:, stage);
    Rate = Rate(Rate > 0);   % drop failed/no-intersection simulations

    MeanRate(stage) = mean(Rate);
    StdRate(stage)  = std(Rate);

    P10(stage) = prctile(Rate, 90);   % 90% probability of exceeding -> conservative
    P50(stage) = prctile(Rate, 50);
    P90(stage) = prctile(Rate, 10);   % 10% probability of exceeding -> optimistic

    CI_lower(stage) = prctile(Rate, 2.5);
    CI_upper(stage) = prctile(Rate, 97.5);

    ProbabilityAbove100(stage) = sum(Rate > 100) / length(Rate);
    ProbabilityAbove80(stage)  = sum(Rate > 80) / length(Rate);
end

MonteCarloSummary = table(PressureCase, MeanRate', StdRate', P10', P50', P90', ...
    CI_lower', CI_upper', ProbabilityAbove100', ProbabilityAbove80', ...
    'VariableNames', {'Case', 'MeanRate', 'StdDev', 'P10', 'P50', 'P90', ...
                       'CI_Lower', 'CI_Upper', 'ProbAbove100', 'ProbAbove80'});
disp(MonteCarloSummary)

% --- Distribution plots ---
figure
PressureLabels = {'Pr=1500', 'Pr=1410', 'Pr=1300', 'Pr=1200'};
for i = 1:nStages
    subplot(2, 2, i)
    histogram(OperatingRate(:, i), 30)
    xlabel('Oil Rate')
    ylabel('Frequency')
    title(['Well A, ', PressureLabels{i}])
    grid on
end

% --- Mean rate & uncertainty bar charts ---
figure
bar(MeanRate)
grid on
xticklabels(PressureLabels)
ylabel('Mean Oil Rate')
title('Expected Production Rate — Well A vs Depletion Stage')

figure
bar(StdRate)
grid on
xticklabels(PressureLabels)
ylabel('Standard Deviation')
title('Production Uncertainty — Well A vs Depletion Stage')

% --- Confidence interval plot ---
figure
errorbar(1:nStages, MeanRate, MeanRate - CI_lower, CI_upper - MeanRate, 'o', 'LineWidth', 2)
grid on
xticks(1:nStages)
xticklabels({'1500', '1410', '1300', '1200'})
xlabel('Reservoir Pressure, psig')
ylabel('Oil Rate')
title('Well A: 95% Confidence Interval vs Depletion Stage')

% --- Probability of meeting targets ---
figure
bar([ProbabilityAbove100; ProbabilityAbove80]')
grid on
legend('Q > 100', 'Q > 80')
xticklabels(PressureLabels)
ylabel('Probability')
title('Well A: Probability of Meeting Production Targets')

% --- Cumulative distribution across stages ---
figure
hold on
for stage = 1:nStages
    Rate = sort(OperatingRate(OperatingRate(:, stage) > 0, stage));
    F = (1:length(Rate)) / length(Rate);
    plot(Rate, F, 'LineWidth', 2)
end
grid on
xlabel('Oil Rate')
ylabel('Cumulative Probability')
legend(PressureLabels)
title('Well A: Monte Carlo Cumulative Distribution Across Depletion Stages')
hold off
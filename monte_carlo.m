%production forecast of a pool
%water production is zero and drive mechanism is depletion drive
OIIP_entire_pool=10*10^6
%current situation:
pavg=1410%psig
oil_produced =669000%bbl

oil = [0,0.031, 0.353, 0.669];
p_avg = [2350,2100, 1720, 1410];

Np = [0.120, 0.229, 0.307, 0.402, 0.471, 0.533, 0.565, 0.602, 0.641, 0.669];
GOR = [209, 208, 214, 220, 242, 240, 255, 298, 353, 365];

figure

% Left Y-axis (Pressure)
yyaxis left
plot(oil, p_avg, '-o','LineWidth',2)
ylabel('Average Pressure')

% Right Y-axis (GOR)
yyaxis right
plot(Np, GOR, '-s','LineWidth',2)
ylabel('GOR')
hold on

% Polynomial regression for GOR vs Np (choose degree, e.g., 2 or 3)
degree = 3; 
p = polyfit(Np, GOR, degree);           % Fit polynomial
Np_smooth = linspace(min(Np), max(Np), 200); % Smooth x values
GOR_fit = polyval(p, Np_smooth);        % Evaluate polynomial
plot(Np_smooth, GOR_fit, 'r--','LineWidth',2) % Plot smooth fit

xlabel('Cumulative Oil Production')
title('Pressure and GOR vs Cumulative Oil')
grid on
legend('Pressure','GOR','GOR Fit','Location','best')

% Data
Px_all = [2350, 1820, 1710, 1420, 2100, 1730,1550, 2100, 1660,1400, 2100, 1770,1420]; % main points only
J_all  = [0.22, 0.19, 0.09, 0.14, 0.06, 0.07,0.05, 0.19, 0.15, 0.12,0.11, 0.09,0.07];

% Normalize by first value of each well manually
jy_all = [1, 0.86, 0.41, 0.64, 1, 1.17,0.83, 1, 0.79,0.63, 1, 0.82,0.64];  % approximate normalized values

% Fit a single straight line through all these points
p = polyfit(Px_all, jy_all, 1);
Px_fit = linspace(min(Px_all), max(Px_all), 200);
jy_fit = polyval(p, Px_fit);

% Plot
figure
hold on
plot(Px_fit, jy_fit, 'k-', 'LineWidth', 2, 'DisplayName', 'Single Trend Line')  % only one line

% Plot markers for all points (optional)
plot(Px_all, jy_all, 'o', 'MarkerSize', 8, 'MarkerFaceColor','r', 'DisplayName','Data Points')

set(gca,'XDir','reverse')  % decreasing X-axis
xlabel('Average static pressure at datum, psig')
ylabel('J / J_i')
title('Production Forecast: Single PI Trend')
grid on
legend('Location','best')
ylim([0 1.5])
hold off


% New x values
pressure = [1500, 1410, 1300, 1200, 1000, 800, 600, 400, 300, 200];

% Use the same fit coefficients 'p' to get y values
jy_pressure = polyval(p, pressure);

% Display results
disp(table(pressure', jy_pressure', 'VariableNames', {'Pressure_psig', 'J_over_Ji'}))
j=jy_pressure

JA=j.*0.22
JB=j.*0.06
JC=j.*0.19
JD=j.*0.11
qa=JA.*pressure
qa_pumping =JA.*(pressure-150)
disp(table(pressure', jy_pressure', ...
    'VariableNames', {'Pressure_psig', 'J_over_Ji'}))


%% ==========================================================
%               MONTE CARLO INITIALIZATION
% ===========================================================

rng(123);                 % For reproducible results

Nsim = 5000;              % Number of simulations

% -----------------------------
% Uncertainty assumptions
% -----------------------------

Pressure_std = 20;        % psi measurement uncertainty
PI_std = 0.10;            % 10% uncertainty in PI
GOR_std = 0.05;           % 5% uncertainty in GOR

% -------------------------------------------------------
% Storage arrays
% -------------------------------------------------------

Pressure_MC = zeros(Nsim,length(pressure));

JA_MC = zeros(Nsim,length(JA));
JB_MC = zeros(Nsim,length(JB));
JC_MC = zeros(Nsim,length(JC));
JD_MC = zeros(Nsim,length(JD));

GOR_MC = zeros(Nsim,length(GOR));

OilRate_MC = zeros(Nsim,length(pressure));

disp('Monte Carlo initialized...')
%% ==========================================================
% Generate Random Reservoir Realizations
% ===========================================================

for sim = 1:Nsim

    % -------------------------------
    % Random reservoir pressure
    % -------------------------------

    Pressure_MC(sim,:) = pressure + Pressure_std*randn(size(pressure));

    % Prevent negative pressure

    Pressure_MC(sim,Pressure_MC(sim,:)<0)=0;

    % -------------------------------
    % Random Productivity Index
    % -------------------------------

    JA_MC(sim,:) = JA .* (1 + PI_std*randn(size(JA)));
    JB_MC(sim,:) = JB .* (1 + PI_std*randn(size(JB)));
    JC_MC(sim,:) = JC .* (1 + PI_std*randn(size(JC)));
    JD_MC(sim,:) = JD .* (1 + PI_std*randn(size(JD)));

    % Prevent negative PI

    JA_MC(sim,JA_MC(sim,:)<0)=0;
    JB_MC(sim,JB_MC(sim,:)<0)=0;
    JC_MC(sim,JC_MC(sim,:)<0)=0;
    JD_MC(sim,JD_MC(sim,:)<0)=0;

    % -------------------------------
    % Random GOR
    % -------------------------------

    GOR_MC(sim,:) = GOR .* (1 + GOR_std*randn(size(GOR)));

    GOR_MC(sim,GOR_MC(sim,:)<0)=0;

    % -------------------------------
    % Oil rate estimate
    % -------------------------------

    max_OilRate_MC(sim,:) = JA_MC(sim,:) .* Pressure_MC(sim,:);

end

disp('Monte Carlo realizations generated.')
%% ==========================================================
% Monte Carlo Results
% ===========================================================

figure

histogram(max_OilRate_MC(:,1),40)

xlabel('Oil Rate')
ylabel('Frequency')

title('Monte Carlo Distribution of Initial Oil Rate')

grid on
figure

histogram(Pressure_MC(:,2),40)

xlabel('Reservoir Pressure (psi)')
ylabel('Frequency')

title('Pressure Uncertainty')

grid on
figure

histogram(JA_MC(:,2),40)

xlabel('Productivity Index')

ylabel('Frequency')

title('PI Uncertainty')

grid on
P10 = prctile(max_OilRate_MC(:,1),90);
P50 = prctile(max_OilRate_MC(:,1),50);
P90 = prctile(max_OilRate_MC(:,1),10);

fprintf('\n');
fprintf('========== Monte Carlo Forecast ==========\n');
fprintf('P10 Oil Rate = %.2f\n',P10);
fprintf('P50 Oil Rate = %.2f\n',P50);
fprintf('P90 Oil Rate = %.2f\n',P90);
fprintf('Mean          = %.2f\n',mean(OilRate_MC(:,1)));
fprintf('Std Dev       = %.2f\n',std(OilRate_MC(:,1)));





% Given
pressure2 = [1500 1410 1300 1200 1000 800 600 400 300 200 0];
Pr = [1500 1410 1300 1200 1000,800,600,400,300,200]                % Reservoir pressure
JA = [JA(1) JA(2) JA(3) JA(4) JA(5) JA(6) JA(7) JA(8) JA(9) JA(10)];  % Your PI values


















%% ============================================
% Monte Carlo IPR Generation (CHANGED: Vogel's Method instead of linear)
% ============================================

Q_MC = cell(Nsim,1);          % Store all IPR curves
PWF_MC = cell(Nsim,1);

Rate_operating = zeros(Nsim,1);

figure
hold on

for sim = 1:Nsim

    JA_random = JA_MC(sim,:);

    for k = 1:length(JA_random)

        Pr_k   = Pressure_MC(sim,k);              % reservoir pressure for this stage/sim
        Qmax_k = JA_random(k) * Pr_k / 1.8;        % Vogel AOF from PI (J-matching relation)

        Pwf_vals = pressure2(k:end);               % same Pwf sweep as before
        Pwf_vals(Pwf_vals > Pr_k) = Pr_k;          % clip so ratio stays in [0,1]

        ratio = Pwf_vals ./ max(Pr_k, eps);        % avoid divide-by-zero
        Q = Qmax_k .* (1 - 0.2*ratio - 0.8*ratio.^2);

        Q(Q<0)=0;

        % Save
        Q_MC{sim,k}=Q;
        PWF_MC{sim,k}=pressure2(k:end);

        % Plot only first 50 simulations
        if sim<=50
            plot(Q, pressure2(k:end),'Color',[0.7 0.7 0.7], 'LineWidth',1)
        end

    end

end
xlabel('Rate')

ylabel('Pressure (psig)')

title('Monte Carlo IPR Curves (Vogel Method)')

grid on


figure

hold on

for sim=1:100

    plot(pressure,JA_MC(sim,:))

end

set(gca,'XDir','reverse')

xlabel('Pressure')

ylabel('PI')

title('Monte Carlo Productivity Index Curves')

grid on


Q_test = [50 100 200];

pwf1 = [990 850 810];
pwf2 = [910 790 710];
pwf3 = [830 730 650];
pwf4 = [780 690 610];

%% ============================================================
% Monte Carlo VLP Curves — Well A, 4 Reservoir Pressure Cases
% (CHANGED: this is ONE well (Well A) at 4 depletion stages,
%  NOT 4 different wells. Loop variable renamed accordingly.)
% ============================================================

VLP_std = 15;          % Pressure uncertainty (psi)

Qfit = linspace(0,250,200);

OperatingRate = zeros(Nsim,4);        % columns = pressure case 1..4 for Well A
OperatingPressure = zeros(Nsim,4);

P10_rate = zeros(1,4);
P50_rate = zeros(1,4);
P90_rate = zeros(1,4);

%% ==========================================
% Monte Carlo Polynomial Fits (Well A, 4 pressure cases)
% ==========================================

for sim = 1:Nsim

    % Random VLP measurements (all Well A test data)

    pwf1_rand = pwf1 + randn(size(pwf1))*VLP_std;
    pwf2_rand = pwf2 + randn(size(pwf2))*VLP_std;
    pwf3_rand = pwf3 + randn(size(pwf3))*VLP_std;
    pwf4_rand = pwf4 + randn(size(pwf4))*VLP_std;

    % Fit curves

    poly{1}=polyfit(Q_test,pwf1_rand,2);
    poly{2}=polyfit(Q_test,pwf2_rand,2);
    poly{3}=polyfit(Q_test,pwf3_rand,2);
    poly{4}=polyfit(Q_test,pwf4_rand,2);


    % CHANGED: loop over pressure CASES for Well A, not different wells
    for caseIdx = 1:4

        Pr_k   = Pressure_MC(sim,caseIdx);            % Well A reservoir pressure at this depletion stage
        Qmax_k = JA_MC(sim,caseIdx) * Pr_k / 1.8;      % Well A Vogel AOF at this stage

        coef = poly{caseIdx};

        Q_grid = linspace(0.001, Qmax_k, 500);
        Pwf_vogel = Pr_k * 0.125 * (-1 + sqrt(81 - 80*(Q_grid/Qmax_k)));
        Pwf_vlp   = polyval(coef, Q_grid);

        diff_curve = Pwf_vogel - Pwf_vlp;
        idx = find(diff(sign(diff_curve)) ~= 0, 1);

        if isempty(idx)
            continue
        end

        Q1 = Q_grid(idx);     Q2 = Q_grid(idx+1);
        d1 = diff_curve(idx); d2 = diff_curve(idx+1);
        q = Q1 - d1*(Q2-Q1)/(d2-d1);
        p_op = polyval(coef, q);

        OperatingRate(sim,caseIdx)=q;
        OperatingPressure(sim,caseIdx)=p_op;

    end

end
%% ==========================================
% Statistics
% ==========================================

for caseIdx=1:4

    rate = OperatingRate(:,caseIdx);

    rate = rate(rate>0);

    P10_rate(caseIdx)=prctile(rate,90);

    P50_rate(caseIdx)=prctile(rate,50);

    P90_rate(caseIdx)=prctile(rate,10);

end

% CHANGED: labels now reflect Well A at 4 reservoir pressure stages,
% not 4 different wells

PressureCase = ["Pr = 1500 psig"; "Pr = 1410 psig"; "Pr = 1300 psig"; "Pr = 1200 psig"];

MonteCarloTable = table(PressureCase,...
    P10_rate',...
    P50_rate',...
    P90_rate',...
    'VariableNames',...
    {'Case','P10','P50','P90'});

disp(MonteCarloTable)
figure

histogram(OperatingRate(:,1),40)

xlabel('Oil Rate')

ylabel('Frequency')

title('Well A Monte Carlo Production Forecast at Pr = 1500 psig')

grid on
figure

PressureLabels = {'Pr=1500','Pr=1410','Pr=1300','Pr=1200'};

for i = 1:4

    subplot(2,2,i)

    histogram(OperatingRate(:,i),30)

    xlabel('Oil Rate')

    ylabel('Frequency')

    title(['Well A, ', PressureLabels{i}])   % CHANGED from "Well ",char('A'+i-1)

    grid on

end




%% ==========================================================
% PART 5 : MONTE CARLO STATISTICS (Well A, 4 pressure cases)
% ===========================================================

MeanRate = zeros(1,4);
StdRate = zeros(1,4);

P5 = zeros(1,4);
P10 = zeros(1,4);
P50 = zeros(1,4);
P90 = zeros(1,4);
P95 = zeros(1,4);

CI_lower = zeros(1,4);
CI_upper = zeros(1,4);

ProbabilityAbove100 = zeros(1,4);
ProbabilityAbove80 = zeros(1,4);

for caseIdx = 1:4

    Rate = OperatingRate(:,caseIdx);

    % Remove failed simulations
    Rate = Rate(Rate>0);

    MeanRate(caseIdx) = mean(Rate);

    StdRate(caseIdx) = std(Rate);

    P5(caseIdx) = prctile(Rate,95);

    P10(caseIdx) = prctile(Rate,90);

    P50(caseIdx) = prctile(Rate,50);

    P90(caseIdx) = prctile(Rate,10);

    P95(caseIdx) = prctile(Rate,5);

    % 95% Confidence Interval

    CI_lower(caseIdx) = prctile(Rate,2.5);

    CI_upper(caseIdx) = prctile(Rate,97.5);

    % Probability of exceeding production targets

    ProbabilityAbove100(caseIdx) = sum(Rate>100)/length(Rate);

    ProbabilityAbove80(caseIdx) = sum(Rate>80)/length(Rate);

end

%% ==========================================================
% Summary Table (Well A, 4 pressure cases)
% ===========================================================

PressureCase = ["Pr = 1500 psig"; "Pr = 1410 psig"; "Pr = 1300 psig"; "Pr = 1200 psig"];

MonteCarloSummary = table( ...
    PressureCase,...
    MeanRate',...
    StdRate',...
    P10',...
    P50',...
    P90',...
    CI_lower',...
    CI_upper',...
    ProbabilityAbove100',...
    ProbabilityAbove80',...
    'VariableNames',...
    {'Case',...
    'MeanRate',...
    'StdDev',...
    'P10',...
    'P50',...
    'P90',...
    'CI_Lower',...
    'CI_Upper',...
    'ProbAbove100',...
    'ProbAbove80'});

disp(MonteCarloSummary)


figure

bar(MeanRate)

grid on

xticklabels({'Pr=1500','Pr=1410','Pr=1300','Pr=1200'})   % CHANGED from Well A/B/C/D

ylabel('Mean Oil Rate')

title('Expected Production Rate — Well A vs Depletion Stage')

figure

bar(StdRate)

grid on

xticklabels({'Pr=1500','Pr=1410','Pr=1300','Pr=1200'})

ylabel('Standard Deviation')

title('Production Uncertainty — Well A vs Depletion Stage')

figure

hold on

errorbar(1:4,...
    MeanRate,...
    MeanRate-CI_lower,...
    CI_upper-MeanRate,...
    'o','LineWidth',2)

grid on

xticks(1:4)

xticklabels({'1500','1410','1300','1200'})

xlabel('Reservoir Pressure, psig')

ylabel('Oil Rate')

title('Well A: 95% Confidence Interval vs Depletion Stage')

hold off

figure

bar([ProbabilityAbove100;
     ProbabilityAbove80]')

grid on

legend('Q > 100','Q > 80')

xticklabels({'Pr=1500','Pr=1410','Pr=1300','Pr=1200'})

ylabel('Probability')

title('Well A: Probability of Meeting Production Targets')

figure

hold on

for caseIdx = 1:4

    Rate = OperatingRate(:,caseIdx);
    Rate = sort(Rate(Rate>0));

    F = (1:length(Rate))/length(Rate);

    plot(Rate,F,'LineWidth',2)

end

grid on

xlabel('Oil Rate')

ylabel('Cumulative Probability')

legend('Pr=1500','Pr=1410','Pr=1300','Pr=1200')

title('Well A: Monte Carlo Cumulative Distribution Across Depletion Stages')



% H+ Diffusion in clinopyroxene 
% Model: Non-isothermal finite slab in an infinite reservoir with non-zero surface concentration

%% Input values 
C0_UL = 275; % Initial H2O concentration (Upper limit)
C0_LL = 192; % Initial H2O concentration (Lower limit)
C1 = 5.5;  % Final H2O concentration (at grain rim)

%h_UL = 579.2116; Input max. half-length of SIMS profile
%h_LL = 460.538; Input min. half-length of SIMS profile

% Grid step size
x = 0:5:(h_UL/2); 

% Cooling rate 
del_T = 1; % cooling rate per hour
del_T_s = del_T/3600; % cooling rate per second

% Monte Carlo settings
nMC = 1000;

% Random inputs from the Monte Carlo values
C0_random = C0_LL + (C0_UL - C0_LL).*rand(nMC,1);
h_random  = h_LL + (h_UL - h_LL).*rand(nMC,1);

% Time increments (seconds)
t_increments = [504, 18540, 33984];

% Diffusion parameters (Bissbort et al. 2022)
D0 = 5.47e-8; 
Ea = 115640;  
R = 8.314; 

%% Temperature setup
T0 = 1100 + 273; % Celsius to Kelvin

% Monte Carlo storage
C_xt_MC = nan(length(x), length(t_increments), nMC);

%% Monte Carlo loop
for imc = 1:nMC

    h_this = h_random(imc);

    % Initial profile (fixed for all time steps)
    C_initial = C0_random(imc) * ones(size(x));

    for it = 1:length(t_increments)

        % Cumulative time
        elapsed_time = sum(t_increments(1:it));

        % Temperature evolution
        T = T0 - elapsed_time * del_T_s;

        % Diffusion coefficient using Arrhenius relation
        D = D0 * exp(-Ea / (R * T));
        D = D * 1e12; % um^2/s

        % New profile for this time step (using the sum of infinite series from j= 0 to 20)
        C_profile = nan(size(x));

        for ix = 1:length(x)
            xx = x(ix);

            if xx > h_this/2
                C_profile(ix) = NaN;
                continue
            end

            sum_terms = 0; 
            for j = 0:20
                n = 2*j + 1;
                term = (1/n) * sin(n*pi*xx/h_this) * ...
                       exp(-((n*pi/h_this)^2) * D * elapsed_time);
                sum_terms = sum_terms + term;
            end

            C_profile(ix) = C1 + (4*(C_initial(ix)-C1)/pi) * sum_terms; 
        end 

        C_xt_MC(:,it,imc) = C_profile;
    end
end

%% Statistics of the Monte Carlo run
C_xt_mean = nanmean(C_xt_MC, 3);
C_xt_min  = nanmin(C_xt_MC, [], 3);
C_xt_max  = nanmax(C_xt_MC, [], 3);

%% Plotting the diffusion profiles
%data = readmatrix('1C_Crystalwise Diffusion Profile.xlsx', 'Sheet','measuredH2O_2a');

x_meas = data(:,4);
C_meas = data(:,5);
C_err  = data(:,6);
group  = data(:,8);

figure;
hold on;

colors = lines(length(t_increments));
t_cumulative = cumsum(t_increments);

h_mean = gobjects(length(t_increments),1);
h_leg  = gobjects(length(t_increments),1);

for pt = 1:length(t_increments)

    % Monte Carlo curves
    for imc = 1:nMC
        plot(x, C_xt_MC(:,pt,imc), ...
            'Color', [colors(pt,:) 0.08]);
    end

    % Envelope
    fill([x fliplr(x)], ...
         [C_xt_min(:,pt)' fliplr(C_xt_max(:,pt)')], ...
         colors(pt,:), ...
         'FaceAlpha', 0.18, ...
         'EdgeColor', 'none');

    % Mean line
    h_mean(pt) = plot(x, C_xt_mean(:,pt), ...
        'LineWidth', 3.5, ...
        'Color', colors(pt,:));

    uistack(h_mean(pt), 'top');

    % Dummy legend line
    h_leg(pt) = plot(nan, nan, ...
        'LineWidth', 3, ...
        'Color', colors(pt,:));
end

%% Overlaying measured cpx H2O data
color_map = [
    0 0 0.5;    
    1 1 0;      
    1 0 0;      
    0.6 1 0.6   
];

unique_groups = unique(group(~isnan(group)));
nGroups = numel(unique_groups);

h_groups = gobjects(nGroups,1);
labels_groups = cell(nGroups,1);

for i = 1:nGroups
    g = unique_groups(i);
    idx = (group == g);

    h_groups(i) = errorbar(x_meas(idx), C_meas(idx), C_err(idx), 'o', ...
        'Color', 'k', ...
        'MarkerFaceColor', color_map(g,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, ...
        'LineStyle', 'none');

    labels_groups{i} = num2str(g);
end

%% Putting axis labels
xlabel('Distance from rim (µm)', 'FontSize', 18, 'Color', 'k');
ylabel('H_{2}O concentration (ppm)', 'FontSize', 18, 'Color', 'k');

title(['Nonisothermal diffusion profile in finite slab (Cooling rate = ', ...
    num2str(del_T), ' °C/hr)'], 'Color', 'k');

%% Legend 1: Diffusion timesteps
labels = arrayfun(@(tt) sprintf('t = %.2f hours', tt/3600), ...
    t_cumulative, 'UniformOutput', false);

lgd1 = legend(h_leg, labels, ...
    'Location','east', 'FontSize', 14, 'TextColor', 'k');
title(lgd1, 'Time', 'Color', 'k');

%% Legend 2: Profile numbers
ax1 = gca;

ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XTick', [], 'YTick', [], ...
           'Box', 'off');

lgd2 = legend(ax2, h_groups, labels_groups, ...
    'Location','northwest', 'FontSize', 14, 'TextColor', 'k');
title(lgd2, 'Profile', 'Color', 'k');

linkaxes([ax1 ax2]);

%% Plot formatting
set(ax1, 'FontSize', 16, ...
         'XColor', 'k', ...
         'YColor', 'k');

ax1.XAxis.FontSize = 14;
ax1.YAxis.FontSize = 14;
ax2.XColor = 'none';
ax2.YColor = 'none';

ylim([0, 280]);
grid off;

saveas(gcf, 'Nonisothermal_DiffusionProfiles.png');

% Parameters
C_source = 1;       % Concentration at water table
D_eff = 0.003;       % Effective diffusion coefficient in cm^2/s
L = 746;              % Depth to surface in cm
z = linspace(0, L, 100); % Discretize the depth

% Time steps (e.g., simulate up to 500 days)
time = linspace(0, 3000*24*3600, 100); % Time in seconds

% Solution matrix
C = zeros(length(z), length(time));

% Calculate concentration over time
for t = 1:length(time)
    C(:, t) = C_source * (1 - erf((L - z) / sqrt(4 * D_eff * time(t))) + erf((L + z) / sqrt(4 * D_eff * time(t))));
end

% Plot the concentration profile over time
figure;
for t = 1:10:length(time)
    plot( C(:, t),z);
    hold on;
end
xlabel('Depth (m)');
ylabel('Concentration');
title('Concentration Profile Over Time');
legend(arrayfun(@(t) sprintf('%.1f days', t/24/3600), time(1:10:end), 'UniformOutput', false));


Time_Eq = L^2/D_eff/(24*3600); 
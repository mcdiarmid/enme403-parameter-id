tutorial_setup

% Initialize forward simulation variables
T0 = mean(T_data(1:2));
dt_ratio = 25;
dt_HR = dt_data/dt_ratio;
t_HR = t_start:dt_HR:t_end;
n_HR = t_end/dt_HR;
phi = zeros(size(t_HR'));

for c = 1:1:3
  for ii = 1:1:n_HR
    if t_impulse(c) <= t_HR(ii)
      phi(ii, c) = 1/dt_HR;
      break
    end
  end
end

% Forward simulation function
Temp_f = @(ts, t0, th, ph) exp(-th(1)*ts) .* (t0+cumtrapz(ts, ...
  exp(th(1)*ts).*((ph*th(2:end))+(th(1)*T_amb))'));

% Initialize
generations = 250;
K = 100;
N = 10;
mut_rate = 0.05;
thetas = 30*rand(K, 5);
max_err = zeros(1, generations);

% Create plots
figure(501);

% K plot
ax_k = subplot(231); hold all;
xlabel(ax_k, 'Generations'); ylabel(ax_k, '\Theta_{opt}(k) [s^{-1}]');
plot(ax_k, [0 generations], [parent(1) parent(1)], '-g');  
xlim(ax_k, [0, generations]); ylim(ax_k, [0 6]);

% J1, J2, J3 plot
ax_j = subplot(232); hold all;
xlabel(ax_j, 'Generations'); ylabel(ax_j, '\Theta_{opt}(J_x) [kJ]');
plot(ax_j, [0 generations], [parent(2) parent(2)]/1000, '-.r'); 
plot(ax_j, [0 generations], [parent(3) parent(3)]/1000, '-.b'); 
plot(ax_j, [0 generations], [parent(4) parent(4)]/1000, '-.m'); 
xlim(ax_j, [0, generations]); ylim(ax_j, [0 600]);

% Residuals plot
ax_ps = subplot(233); hold all;
xlabel(ax_ps, 'Generations'); ylabel(ax_ps, '\Psi_{opt}');
xlim(ax_k, [0, generations]);

% Forward simulation plot with theta opt
ax_sim = subplot(212); hold all;
xlabel(ax_sim, 'Time [s]'); ylabel(ax_sim, 'Temperature [^oC]'); 
ylim(ax_sim, [20 45]);

% Loop for each generation
for gen = 1:generations
  
  % Iterate through parameter sets
  for ii = 1:K
    % Determine psi of theta(:, i)
    theta = thetas(ii, 2:end)';
    T_fwd = Temp_f(t_HR, T0, theta, phi);
    thetas(ii, 1) = norm(T_fwd(1:dt_ratio:end)' - T_data');
  end
  
  % Sort thetas by the lowest (best) residuals
  thetas = sortrows(thetas);
  
  % Check for NaNs
  for ii = 1:K
    if isnan(thetas(ii, 1)) || isnan(thetas(ii, 1))
      % Create an iterator of length K for only "good" organisms
      ind=1; jj=1; iterator = zeros(1, K);
      while jj <= 100
        iterator(jj) = ind;
        jj = jj + 1;
        if mod(jj, ceil(K/ii)) == 0
         ind=ind+1;
        end
      end
      
      % Appropriately remove NaN and Inf sets, replace with better clones
      thetas(:, :) = thetas(iterator, :);
      break
    end
  end
  
  % Discard the worst and clone multiple copies of the top ranking sets
  thetas(:, :) = thetas(ceil(((1:K)/20).^2), :);

  % Mutate, keep one copy of  
  thetas(2:end, 2:end) = thetas(2:end, 2:end) + ...
    mut_rate*[thetas(2:end,1) thetas(2:end,1) thetas(2:end,1) thetas(2:end,1)].*randn(K-1,4);
  
  % Plot of theta opt
  if mod(gen, 10) == 1
    T_fwd = Temp_f(t_HR, T0, thetas(1, 2:end)', phi);
    plot(ax_sim, t_HR, T_fwd, 'r');
  end
  
  % Plot optimal values from this generation
  ps = thetas (1,1);
  
  k = thetas(1,2); J = thetas(1,3:end)*V*Cp;
  fprintf('k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, J(1), J(2), J(3));
  
  plot(ax_k, gen, k, '.g');
  plot(ax_j, gen, J(1)/1000, '.r', gen, J(2)/1000, '.b', gen, J(3)/1000, '.m');
  plot(ax_ps, gen, ps, '.k');
end


plot(ax_sim, t_HR, T_fwd, 'k');
plot(ax_sim, t_meas, T_data, '+b');
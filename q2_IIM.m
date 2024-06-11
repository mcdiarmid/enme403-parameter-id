% Set-up
tutorial_setup

n_loops = 100;
dt_ratio = 10;
dt_HR = dt_data/dt_ratio;  % High resolution time delta
t_HR = t_start:dt_HR:t_end;  % High resolution time array
[~, N_HR] = size(t_HR);

% Matrix initialization
phi = zeros(N_HR, cols-1);
psi = zeros(1, n_loops);

% Populate A and b matrix
T0 = mean(T_data(1:2));

for c = 1:1:3
  for ii = 1:1:N_HR
    if t_impulse(c) <= t_HR(ii)
      phi(ii, c) = 1/dt_HR;
      A(ceil(ii/dt_ratio)+1:end, c+1) = 1;
      break
    end
  end
end


% Simulate and iterate
set(0,'defaultfigureposition',[60 60 420,330])
figure(1)
set(gcf,'position',[60 60 680 330])
h1=axes('position',[.085 .1 .40 .85]);hold all
h2=axes('position',[.585 .1 .40 .85]);hold all

figure(2)
h3=axes;hold all
yyaxis left; hold all;
plot(h3, [0 n_loops], [parent(2) parent(2)]/1000, 'b-');
plot(h3, [0 n_loops], [parent(3) parent(3)]/1000, 'g-');
plot(h3, [0 n_loops], [parent(4) parent(4)]/1000, 'r-');
yyaxis right; hold all;
plot(h3, [0 n_loops], [parent(1) parent(1)], 'k-');
plot(h1, t_meas, T_data, '+k');

% Populate b matrix and make initial guess
T0 = mean(T_data(1:2));
T_fwd = zeros(size(t_HR)); % t_HR = 0:0.025:10 [1x401]
b = T_data' - T0;

for ii = 1:n_loops  % 1:100
  % Update A column 1 (only one that changes), calculate theta
  iTg = cumtrapz(t_HR, -(T_fwd - T_amb));  
  A(:, 1) = iTg(1:dt_ratio:end)';  % dt_ratio = 10 -> iTg(1:10:401)' [41x1]
  theta = A\b;

  % Calculate parameter values and print
  k = theta(1);
  J1 = theta(2)*V*Cp; J2 = theta(3)*V*Cp; J3 = theta(4)*V*Cp;
  fprintf('k=%.2f, J1=%.0f, J2=%.0f, J3=%.0f\n', k, J1, J2, J3);

  % Forward simulate and store residuals
  % Analytical
  T_fwd = exp(-k*t_HR).*(T0+cumtrapz(t_HR, exp(k*t_HR).* ...
    ((phi * theta(2:end)) + (k*T_amb))'));
  
  % Picard: Error stepping, doesnt work
%   T_fwd = T0 + cumtrapz(t_HR, -k*(T_fwd-T_amb) + (phi*theta(2:end))');
%   % Works but doesnt make sense
%   for jj = 2:N_HR  % 2:401
%     T_fwd(1:jj) = T0 + cumtrapz(t_HR(1:jj), -k*(T_fwd(1:jj)-T_amb) + (phi(1:jj, :)*theta(2:end))');
%   end

  % Time-stepping: Euler 1st order
%   T_fwd(1) = T0;
%   for jj = 2:N_HR  % 2:401
%     T_fwd(jj) = T_fwd(jj-1) + (-k*(T_fwd(jj-1)-T_amb) + (phi(jj-1,:)*theta(2:end))')*dt_HR;
%   end

  % Store residuals
  psi(ii) = norm(T_fwd(1:dt_ratio:end)-T_data);

  % Plot forward simulation every 5th iteration
  if mod(ii, 5) == 0
    plot(h1, t_HR, T_fwd, 'r:');ylim(h1, [20, 45]); 
    ylabel(h1, 'Temp (^oC)'); xlabel(h1, 'Time (s)');
  end
  
  % Plot Residuals
  plot(h2, ii, psi(ii), 'ob','markerfacecolor','b');
  ylabel(h2, '\bfresiduals \it\psi'); xlabel(h2, 'Iterations'); 
  
  % Plot paramater convergence
  yyaxis left;  % J1, J2, J3 on left axis
  plot(h3, ii, J1/1000, '.b', ii, J2/1000, '.g', ii, J3/1000, '.r');
  ylabel(h3, 'J_x (kJ)'); xlabel(h3, 'Iterations'); ylim(h3, [0, 600]);
  
  yyaxis right; % k on right axis
  plot(h3, ii, k, '.k');
  ylabel(h3, 'k (s^-^1)'); xlabel(h3, 'Iterations'); ylim(h3, [0, 5]);
end

plot(h1, t_HR, T_fwd, 'b.');ylim(h1, [20, 45]); 
ylabel(h1, 'Temp (^oC)'); xlabel(h1, 'Time (s)');
